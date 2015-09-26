/*
    Copyright (C) 2010  Francesc Alted
    http://blosc.pytables.org
    License: MIT (see LICENSE.txt)

    Filter program that allows the use of the Blosc filter in HDF5.

    This is based on the LZF filter interface (http://h5py.alfven.org)
    by Andrew Collette.

*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include "hdf5.h"
#include "blosc_filter.h"

#if H5Epush_vers == 2
/* 1.8.x */
#if defined(__GNUC__)
#define PUSH_ERR(func, minor, str, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, str, ##__VA_ARGS__)
#elif defined(_MSC_VER)
#define PUSH_ERR(func, minor, str, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, str, __VA_ARGS__)
#else
/* This version is portable but it's better to use compiler-supported
   approaches for handling the trailing comma issue when possible. */
#define PUSH_ERR(func, minor, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, __VA_ARGS__)
#endif	/* defined(__GNUC__) */
#else
/* 1.6.x */
#define PUSH_ERR(func, minor, str) H5Epush(__FILE__, func, __LINE__, H5E_PLINE, minor, str)
#endif

#if H5Pget_filter_by_id_vers == 2
/* 1.8.x */
#define GET_FILTER(a,b,c,d,e,f,g) H5Pget_filter_by_id(a,b,c,d,e,f,g,NULL)
#else
/* 1.6.x */
#define GET_FILTER H5Pget_filter_by_id
#endif

#if H5Z_class_t_vers == 2
/* 1.8.x where x >= 3 */
#define H5Z_16API 0
#else
/* 1.6.x and 1.8.x with x < 3*/
#define H5Z_16API 1
#endif

/* Function prototypes */

size_t blosc_filter(unsigned flags, size_t cd_nelmts,
                    const unsigned cd_values[], size_t nbytes,
                    size_t *buf_size, void **buf);

herr_t blosc_set_local(hid_t dcpl, hid_t type, hid_t space);


/* Variables for managing thread-specific blosc contexts.
   Using a context specific to each calling thread avoids the need
   to synchronize access to blosc calls with a global lock. */

/**
 * Key for getting/setting the thread-local blosc context.
 */
static pthread_key_t blosc_ctx_tls_key;
/**
 * Used to ensure \link #blosc_ctx_tls_key is only initialized once.
 */
static pthread_once_t blosc_ctx_tls_key_once = PTHREAD_ONCE_INIT;

/* TEMP: Provide this typedef from blosc.h */
typedef int blosc_context;
typedef blosc_context* blosc_context_handle;

/**
 * Destroys a thread-local blosc context handle when a thread is exiting.
 */
static void
destroy_blosc_ctx(void* const ptr)
{
    /* 'ptr' is actually the blosc context handle. */
    blosc_context_handle blosc_ctx_handle = (blosc_context_handle)ptr;

    /* TODO: Re-implement this correctly using functions from blosc.h.
       The code below is just a placeholder for testing purposes. */
    free(blosc_ctx_handle);
}


static void
create_blosc_ctx_tls_key()
{
    /* TODO: Check the return value of pthread_key_create, and if it's non-zero,
       report the error in some way (perhaps via the HDF5 macros). */
    (void)pthread_key_create(&blosc_ctx_tls_key, destroy_blosc_ctx);
}

/**
 * Get the blosc context handle from the thread-local storage slot.
 * If the context hasn't yet been created for this thread, it is created
 * before returning.
 */
static int
get_blosc_ctx_tls(blosc_context_handle* blosc_ctx_handle)
{
    int result;
    blosc_context_handle handle;
    if ((handle = (blosc_context_handle)pthread_getspecific(blosc_ctx_tls_key)) == NULL) {
        /* The blosc context hasn't been created for this thread, so create it. */
        handle = (blosc_context_handle)malloc(sizeof(blosc_context)); /* TODO: Create this via some blosc call.
                                                                          This malloc is just a placeholder! */

        /* Store the context handle to the thread-local storage slot. */
        if ((result = pthread_setspecific(blosc_ctx_tls_key, (void*)handle)) != 0) {
            *blosc_ctx_handle = NULL;
            return result;
        }
    }

    /* Store the context handle before returning. */
    *blosc_ctx_handle = handle;
    return 0;
}

/* Register the filter, passing on the HDF5 return value */
int register_blosc(char **version, char **date){

    int retval;

#if H5Z_16API
    H5Z_class_t filter_class = {
        (H5Z_filter_t)(FILTER_BLOSC),
        "blosc",
        NULL,
        (H5Z_set_local_func_t)(blosc_set_local),
        (H5Z_func_t)(blosc_filter)
    };
#else
    H5Z_class_t filter_class = {
        H5Z_CLASS_T_VERS,
        (H5Z_filter_t)(FILTER_BLOSC),
        1, 1,
        "blosc",
        NULL,
        (H5Z_set_local_func_t)(blosc_set_local),
        (H5Z_func_t)(blosc_filter)
    };
#endif

    retval = H5Zregister(&filter_class);
    if(retval<0){
        PUSH_ERR("register_blosc", H5E_CANTREGISTER, "Can't register Blosc filter");
    }
    *version = strdup(BLOSC_VERSION_STRING);
    *date = strdup(BLOSC_VERSION_DATE);

    /* Initialize the thread-local storage key for blosc contexts. */
    (void) pthread_once(&blosc_ctx_tls_key_once, create_blosc_ctx_tls_key);

    return 1; /* lib is available */
}

/*  Filter setup.  Records the following inside the DCPL:

    1. If version information is not present, set slots 0 and 1 to the filter
       revision and Blosc version, respectively.

    2. Compute the type size in bytes and store it in slot 2.

    3. Compute the chunk size in bytes and store it in slot 3.
*/
herr_t blosc_set_local(hid_t dcpl, hid_t type, hid_t space){

    int ndims;
    int i;
    herr_t r;

    unsigned int typesize, basetypesize;
    unsigned int bufsize;
    hsize_t chunkdims[32];
    unsigned int flags;
    size_t nelements = 8;
    unsigned int values[] = {0,0,0,0,0,0,0,0};
    hid_t super_type;
    H5T_class_t classt;

    r = GET_FILTER(dcpl, FILTER_BLOSC, &flags, &nelements, values, 0, NULL);
    if(r<0) return -1;

    if(nelements < 4) nelements = 4;  /* First 4 slots reserved. */

    /* Set Blosc info in first two slots */
    values[0] = FILTER_BLOSC_VERSION;
    values[1] = BLOSC_VERSION_FORMAT;

    ndims = H5Pget_chunk(dcpl, 32, chunkdims);
    if(ndims<0) return -1;
    if(ndims>32){
        PUSH_ERR("blosc_set_local", H5E_CALLBACK, "Chunk rank exceeds limit");
        return -1;
    }

    typesize = H5Tget_size(type);
    if (typesize==0) return -1;
    /* Get the size of the base type, even for ARRAY types */
    classt = H5Tget_class(type);
    if (classt == H5T_ARRAY) {
      /* Get the array base component */
      super_type = H5Tget_super(type);
      basetypesize = H5Tget_size(super_type);
      /* Release resources */
      H5Tclose(super_type);
    }
    else {
      basetypesize = typesize;
    }

    /* Limit large typesizes (they are pretty inneficient to shuffle
       and, in addition, Blosc does not handle typesizes larger than
       blocksizes). */
    if (basetypesize > BLOSC_MAX_TYPESIZE) basetypesize = 1;
    values[2] = basetypesize;

    /* Get the size of the chunk */
    bufsize = typesize;
    for (i=0; i<ndims; i++) {
        bufsize *= chunkdims[i];
    }
    values[3] = bufsize;

#ifdef BLOSC_DEBUG
    fprintf(stderr, "Blosc: Computed buffer size %d\n", bufsize);
#endif

    r = H5Pmodify_filter(dcpl, FILTER_BLOSC, flags, nelements, values);
    if(r<0) return -1;

    return 1;
}


/* The filter function */
size_t blosc_filter(unsigned flags, size_t cd_nelmts,
                    const unsigned cd_values[], size_t nbytes,
                    size_t *buf_size, void **buf){

    void* outbuf = NULL;
    int status = 0;                /* Return code from Blosc routines */
    size_t typesize;
    size_t outbuf_size;
    int clevel = 5;                /* Compression level default */
    int doshuffle = 1;             /* Shuffle default */
    int compcode;                  /* Blosc compressor */
    int code;
    char *compname = "blosclz";    /* The compressor by default */
    char *complist;
    char errmsg[256];

    /* Filter params that are always set */
    typesize = cd_values[2];      /* The datatype size */
    outbuf_size = cd_values[3];   /* Precomputed buffer guess */
    /* Optional params */
    if (cd_nelmts >= 5) {
        clevel = cd_values[4];        /* The compression level */
    }
    if (cd_nelmts >= 6) {
        doshuffle = cd_values[5];     /* Shuffle? */
    }
    if (cd_nelmts >= 7) {
        compcode = cd_values[6];     /* The Blosc compressor used */
	/* Check that we actually have support for the compressor code */
        complist = blosc_list_compressors();
	code = blosc_compcode_to_compname(compcode, &compname);
	if (code == -1) {
#if H5Epush_vers == 2
            PUSH_ERR("blosc_filter", H5E_CALLBACK,
                     "this Blosc library does not have support for "
                     "the '%s' compressor, but only for: %s",
                     compname, complist);
#else
	    sprintf(errmsg, "this Blosc library does not have support for "
                    "the '%s' compressor, but only for: %s",
		    compname, complist);
            PUSH_ERR("blosc_filter", H5E_CALLBACK, errmsg);
            goto failed;
#endif
	}
    }

    /* We're compressing */
    if(!(flags & H5Z_FLAG_REVERSE)){

#ifdef BLOSC_DEBUG
        fprintf(stderr, "Blosc: Compress %zd chunk w/buffer %zd\n",
		nbytes, outbuf_size);
#endif

        /* Allocate an output buffer exactly as long as the input data; if
           the result is larger, we simply return 0.  The filter is flagged
           as optional, so HDF5 marks the chunk as uncompressed and
           proceeds.
        */

        outbuf_size = (*buf_size);
        outbuf = malloc(outbuf_size);

        if(outbuf == NULL){
            PUSH_ERR("blosc_filter", H5E_CALLBACK,
                     "Can't allocate compression buffer");
            goto failed;
        }

#if ( (BLOSC_VERSION_MAJOR <= 1) && (BLOSC_VERSION_MINOR < 5) )
        status = blosc_compress(clevel, doshuffle, typesize, nbytes,
                                *buf, outbuf, nbytes);
#else
        /* Starting from Blosc 1.5 on, there is not an internal global
	   lock anymore, so do not try to run in multithreading mode
	   so as to not interfering with other possible threads
	   launched by the main Python application */
        status = blosc_compress_ctx(clevel, doshuffle, typesize, nbytes,
                                    *buf, outbuf, nbytes, compname, 0, 1);
#endif
        if (status < 0) {
          PUSH_ERR("blosc_filter", H5E_CALLBACK, "Blosc compression error");
          goto failed;
        }

    /* We're decompressing */
    } else {
        /* declare dummy variables */
        size_t cbytes, blocksize;

#ifdef BLOSC_DEBUG
        fprintf(stderr, "Blosc: Decompress %zd chunk w/buffer %zd\n", nbytes, outbuf_size);
#endif

        free(outbuf);

        /* Extract the exact outbuf_size from the buffer header.
         *
         * NOTE: the guess value got from "cd_values" corresponds to the
         * uncompressed chunk size but it should not be used in a general
         * cases since other filters in the pipeline can modify the buffere
         *  size.
         */
        blosc_cbuffer_sizes(*buf, &outbuf_size, &cbytes, &blocksize);

        outbuf = malloc(outbuf_size);

        if(outbuf == NULL){
          PUSH_ERR("blosc_filter", H5E_CALLBACK, "Can't allocate decompression buffer");
          goto failed;
        }

#if ( (BLOSC_VERSION_MAJOR <= 1) && (BLOSC_VERSION_MINOR < 5) )
	status = blosc_decompress(*buf, outbuf, outbuf_size);
#else
        /* Starting from Blosc 1.5 on, there is not an internal global
	   lock anymore, so do not try to run in multithreading mode
	   so as to not interfering with other possible threads
	   launched by the main Python application */
        status = blosc_decompress_ctx(*buf, outbuf, outbuf_size, 1);
#endif

        if(status <= 0){    /* decompression failed */
          PUSH_ERR("blosc_filter", H5E_CALLBACK, "Blosc decompression error");
          goto failed;
        } /* if !status */

    } /* compressing vs decompressing */

    if(status != 0){
        free(*buf);
        *buf = outbuf;
        *buf_size = outbuf_size;
        return status;  /* Size of compressed/decompressed data */
    }

 failed:
    free(outbuf);
    return 0;

} /* End filter function */
