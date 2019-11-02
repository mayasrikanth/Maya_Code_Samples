/*! \file
 * Declarations for a simple memory allocator.  The allocator manages a small
 * pool of memory, provides memory chunks on request, and reintegrates freed
 * memory back into the pool.
 *
 * Adapted from Andre DeHon's CS24 2004, 2006 material.
 * Copyright (C) California Institute of Technology, 2004-2009.
 * All rights reserved.
 */


 struct header;
 struct footer;


/*! Specifies the size of the memory pool the allocator has to work with. */
extern int MEMORY_SIZE;

/* Specifies size of header + size of footer */
extern int HEADER_FOOTER;

/* Initializes allocator state, and memory pool state too. */
void init_myalloc();


/* Attempt to allocate a chunk of memory of "size" bytes. */
unsigned char * myalloc(int size);


/* Free a previously allocated pointer. */
void myfree(unsigned char *oldptr);


/* Clean up the allocator and memory pool state. */
void close_myalloc();

/* Sanity_check function to check total space of memory pool. */
int sanity_check();

 /* Foward coalescing function. */
void forward_coalesce(unsigned char* oldptr);

/* Backward coalescing function. */
void backward_coalesce(unsigned char *oldptr);
