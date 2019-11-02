/*! \file
 * Implementation of a simple memory allocator.  The allocator manages a small
 * pool of memory, provides memory chunks on request, and reintegrates freed
 * memory back into the pool.
 *
 * Adapted from Andre DeHon's CS24 2004, 2006 material.
 * Copyright (C) California Institute of Technology, 2004-2010.
 * All rights reserved.
 * Author: Maya Srikanth
 *         Caltech 2021
 */

#include <stdio.h>
#include <stdlib.h>

#include "myalloc.h"


/*
 * These variables are used to specify the size and address of the memory pool
 * that the simple allocator works against.  The memory pool is allocated within
 * init_myalloc(), and then myalloc() and free() work against this pool of
 * memory that mem points to.
 */
int MEMORY_SIZE;
int HEADER_FOOTER;
unsigned char *mem;

/*
 * The header of each memory block is a 32-bit integer indicating size
 * of the payload (excluding size of header and footer size). A negative integer
 * indicates an occupied space, while a positive integer indicates
* a free space.
 */
typedef struct header
{
  int block_size;   /* Indicates size of following memory block */
} header;

/*
 * The footer of each memory block is a 32-bit integer indicating size
 * of the payload (excluding size of header and footer size). A negative
 * integer indicates an occupied space, while a positive integer indicates
 * a free space.
 */

typedef struct footer
{
  int block_size;
}footer;

/* freeptr is updated during allocation to traverse memory pool. */
static unsigned char *freeptr;


/*!
 * This function initializes both the allocator state, and the memory pool.  It
 * must be called before myalloc() or myfree() will work at all.
 *
 * Note that we allocate the entire memory pool using malloc().  This is so we
 * can create different memory-pool sizes for testing.  Obviously, in a real
 * allocator, this memory pool would either be a fixed memory region, or the
 * allocator would request a memory region from the operating system (see the
 * C standard function sbrk(), for example).
 */
void init_myalloc() {

    /*
     * Allocate the entire memory pool, from which our simple allocator will
     * serve allocation requests.
     */
    mem = (unsigned char *) malloc(MEMORY_SIZE);
    if (mem == 0) {
        fprintf(stderr,
                "init_myalloc: could not get %d bytes from the system\n",
		MEMORY_SIZE);
        abort();
    }
    /*
     * Initializing state of memory pool by assigning each block a header
     * and footer.
      */
    header *initial_header = (header *)mem;
    initial_header->block_size = MEMORY_SIZE - sizeof(header) - sizeof(footer);
    footer *initial_footer  = (footer *)(mem + MEMORY_SIZE - sizeof(footer));
    initial_footer->block_size = MEMORY_SIZE - sizeof(header) - sizeof(footer);

    freeptr = mem;
}


/*
 * myalloc(): traverses the memory pool implicitly and selects a block of closest
 * size to the allocation request, then devoting the free space to the
 * allocation request. Depending on the size of the chosen block, the
 * allocator will split the block or leave it whole.
 *
 * The function implements a best fit strategy. This minimizes fragmentation
 * by ensuring that the closest free block is dedicated to the allocation
 * request. This leaves the bigger free blocks unfragmented, an important
 * strategic advantage in memory allocation. If a big allocation request comes,
 * the request is more likely to be satisfied if fragmentation is minimized.
 * A disadvantage of best fit strategy is linear time...the whole memory pool
 * must be traversed in order to pinpoint the most apt free space, whereas
 * other strategies may allow for constant time selection of an appropriate
 * block.
 * Time complexity: Linear time
 */
unsigned char *myalloc(int size) {


     freeptr = mem;
     header *fit;   /* Best-fitting block. */
     unsigned char *resultptr;
     int left_over = MEMORY_SIZE - size;
     int allocated = 0;
     HEADER_FOOTER = sizeof(header) + sizeof(footer);

     /* Choosing the block that minimizes left_over. */
     while(freeptr < mem + MEMORY_SIZE) {

       if(((header *)freeptr)->block_size >= size)
       {
         if(((header *)freeptr)->block_size - size < left_over) {
           left_over = ((header *)freeptr)->block_size - size;
           fit = (header *)freeptr;
           allocated = 1;
       }
     }
     freeptr = (unsigned char *) ((unsigned char*)freeptr + HEADER_FOOTER
      + abs(((header *)freeptr)->block_size));
   }

   /* Not enough memory was available for allocation request. */
   if(allocated == 0) {
         fprintf(stderr, "myalloc: cannot service request of size %d with"
                 " %lx bytes allocated\n", size, (freeptr - mem));
         return (unsigned char *) 0;
   }
    /*
     * If block_size exactly fits the requested size, no block splitting.
     */
     resultptr = (unsigned char *)fit;
     footer* fit_footer = (footer *)((unsigned char*)fit + sizeof(header)
       + fit->block_size);

    if(fit->block_size < size + HEADER_FOOTER) {
        /* If block is free space, begin allocation process. */

        /* Updating header value to read 'allocated'. */
        fit->block_size = fit->block_size * -1;

        /* Updating footer value to read 'allocated'. */
        fit_footer->block_size = fit->block_size * -1;
        resultptr = resultptr + sizeof(header);

        return resultptr;
    }

  /* Save value of current block_size in 'old_block'.   */
    int old_block = fit->block_size;

  /* Determine amount of free space remaining in block after allocation. */
    int second_segment = old_block - (size + HEADER_FOOTER);

  /* Update block_size value of header to reflect requested size. */
    fit->block_size = size * -1;


  /* Create footer for block of requested size. */
    footer *new_footer = (footer *) ((unsigned char*)fit + size +
      sizeof(header));
    new_footer->block_size = size * -1;

  /* Create header for remaining space in block. */
    header *new_header = (header *) ((unsigned char*)new_footer
      + sizeof(footer));
    new_header->block_size = second_segment;


  /* Update footer of old_block to reflect second_segment size. */
    fit_footer->block_size = second_segment;

   resultptr = resultptr + sizeof(header);
    //header = resultptr - HEADER_SIZE
    //resultptr = resultptr + sizeof(header);
    return resultptr;

}

/* Simple sanity check for memory pool size.     */
 int sanity_check() {

   unsigned char *ptr = mem;
   int total_space = 0;

   while(ptr < mem + MEMORY_SIZE)
   {
     total_space += abs(((header *)ptr)->block_size);
     ptr = (unsigned char*)((unsigned char*)ptr + HEADER_FOOTER);
     ptr = ptr + ptr->block_size;
     total_space += HEADER_FOOTER;
   }
   return total_space;
 }



/*
 * myfree: Calls forward_coalesce and backward_coalesce, freeing a
 * previously allocated pointer and catalyzing forward or backward coalescing
 * when permissible. Time complexity of deallocation is constant,
 * as it involves altering the sign of block_size to reflect state of
 * memory block.
 *
 */
void myfree(unsigned char *oldptr) {

    /* Marking block as free. */
    header *ptr = (header *)((unsigned char*) oldptr - sizeof(header));
    ptr->block_size = -1 * ptr->block_size;
    forward_coalesce(oldptr);
    backward_coalesce(oldptr);

}

/*
 * forward_coalesce: forward-coallesces a newly deallocated block with a
 * contiguous free block following it in memory. Ensures that the next
 * block exists within bounds of the memory pool before coalescing.
 * Time complexity of alloation is constant time, as a pointer to the
 * current block simply incremented an appropriate amount of bytes to
 * access and update the tags on the following block.
 *
 */
void forward_coalesce(unsigned char* oldptr) {

    header *ptr = (header *)((unsigned char*) oldptr - sizeof(header));

    /* Creating foot to point to footer of deallocated block. */
    footer *foot = (footer *) ((unsigned char *)ptr + ptr->block_size) + 1;
    foot->block_size = ptr->block_size;

    /* Ensure next block exists before coalescing with next block.*/
    if((unsigned char *) ptr + ptr->block_size + HEADER_FOOTER + 4
      < mem + MEMORY_SIZE) {

    /* Check to see if following block is free. */
      header *next_block = (header *)((unsigned char *)ptr + ptr->block_size
        + 2 * sizeof(header));
      footer *next_block_footer;

      if(next_block->block_size >= 0) {
        /* Updating header block_size of oldptr block. */
        ptr->block_size = ptr->block_size + next_block->block_size;
        ptr->block_size = ptr->block_size + 2 * sizeof(header);

        /* Storing new_block_size in variable. */
        int new_block_size = ptr->block_size;

        next_block_footer = (footer *)((unsigned char*)next_block
          + next_block->block_size + sizeof(header));

        /* Updating footer block_size of next_block. */
        next_block_footer->block_size = new_block_size;

      }
    }
}

/*
 * backward_coalesce: backward-coallesces a newly deallocated block with a
 * contiguous free block. Ensures that the previous block exists within bounds
 * of the memory pool before coalescing. Time complexity of alloation
 * is constant time, as a pointer to the current block simply decremented
 * an appropriate amount of bytes to access and update the tags on the
 * previous block.
 *
 */

  void backward_coalesce(unsigned char *oldptr) {

    header *ptr = (header *)((unsigned char*) oldptr - sizeof(header));

    /* Ensure previous block exists before coalescing with previous block. */
    if(oldptr - HEADER_FOOTER > mem) {


      /* Ensure previous block is free. */
      footer *prev_block = (footer *)((unsigned char *)ptr - sizeof(footer));
      header *prev_block_header;
      footer *old_ptr_foot;

      if(prev_block->block_size >= 0) {

        prev_block_header = (header *)((unsigned char *)ptr
          - 2 * sizeof(footer) - prev_block->block_size);
        /* Updating header of previous block. */
        prev_block_header->block_size = ptr->block_size +
          prev_block->block_size + 2 * sizeof(header);
        /* Updating footer of oldptr block. */

        old_ptr_foot = (footer *)((unsigned char *)ptr +
          ptr->block_size + sizeof(header));
        old_ptr_foot->block_size = ptr->block_size +
          prev_block->block_size + 2 * sizeof(header);
       }
     }
   }


/*!
 * Clean up the allocator state.
 * All this really has to do is free the user memory pool. This function mostly
 * ensures that the test program doesn't leak memory, so it's easy to check
 * if the allocator does.
 */
void close_myalloc() {
    free(mem);
}
