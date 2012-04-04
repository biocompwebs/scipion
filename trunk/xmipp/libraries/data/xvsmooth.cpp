/*
 * xvsmooth.c - smoothing/color dither routines for XV
 *
 *  Author:    John Bradley, University of Pennsylvania
 *                (bradley@cis.upenn.edu)
 *
 *  Contains:
 *            byte *SmoothResize(src8, swide, shigh, dwide, dhigh,
 *                               rmap, gmap, bmap, rdmap, gdmap, bdmap, maplen)
 *            byte *Smooth24(pic824, swide, shigh, dwide, dhigh,
 *                               rmap, gmap, bmap)
 *            byte *DoColorDither(picSrc, pic8, w, h, rmap,gmap,bmap,
 *                                rdisp, gdisp, bdisp, maplen)


 */

/* Copyright Notice
 * ================
 * Copyright 1989, 1990, 1991, 1992, 1993 by John Bradley
 *
 * Permission to use, copy, and distribute XV in its entirety, for
 * non-commercial purposes, is hereby granted without fee, provided that
 * this license information and copyright notice appear in all copies.
 *
 * Note that distributing XV 'bundled' in with ANY product is considered
 * to be a 'commercial purpose'.
 *
 * Also note that any copies of XV that are distributed MUST be built
 * and/or configured to be in their 'unregistered copy' mode, so that it
 * is made obvious to the user that XV is shareware, and that they should
 * consider donating, or at least reading this License Info.
 *
 * The software may be modified for your own purposes, but modified
 * versions may NOT be distributed without prior consent of the author.
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the author be held liable for any damages
 * arising from the use of this software.
 *
 * If you would like to do something with XV that this copyright
 * prohibits (such as distributing it with a commercial product,
 * using portions of the source in some other program, etc.), please
 * contact the author (preferably via email).  Arrangements can
 * probably be worked out.
 *
 * XV is shareware for PERSONAL USE only.  You may use XV for your own
 * amusement, and if you find it nifty, useful, generally cool, or of
 * some value to you, your non-deductable donation would be greatly
 * appreciated.  $25 is the suggested donation, though, of course,
 * larger donations are quite welcome.  Folks who donate $25 or more
 * can receive a Real Nice bound copy of the XV manual for no extra
 * charge.
 *
 * Commercial, government, and institutional users MUST register their
 * copies of XV, for the exceedingly REASONABLE price of just $25 per
 * workstation/X terminal.  Site licenses are available for those who
 * wish to run XV on a large number of machines.  Contact the author
 * for more details.
 *
 * The author may be contacted via:
 *    US Mail:  John Bradley
 *              1053 Floyd Terrace
 *              Bryn Mawr, PA  19010
 *
 *    Phone:    (215) 898-8813
 *    EMail:    bradley@cis.upenn.edu
 */


/* CO: --------------------------------------------------------------------- */
/* RANGE forces a to be in the range b..c (inclusive) */
#define RANGE(a,b,c) { if (a < b) a = b;  if (a > c) a = c; }
typedef unsigned char byte;

/* CO: --------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "xmipp_error.h"

byte *Smooth(byte *picSrc8, int swide, int shigh, int dwide, int dhigh);
void DoColorDither(byte *picSmooth, byte *&picDithered, int w, int h);

#ifdef __STDC__
int SmoothX(byte *, byte *, int, int, int, int);
int SmoothY(byte *, byte *, int, int, int, int);
int SmoothXY(byte *, byte *, int, int, int, int);
#else
int  SmoothX(), SmoothY(), SmoothXY();
#endif

#define xvbzero(s,l) memset(s, 0, l)
/***************************************************/
byte *SmoothResize(byte *picSrc8, int swide, int shigh,
                   int dwide, int dhigh)
{
    /* generic interface to Smooth and ColorDither code.
       given an 8-bit-per, swide * shigh image with colormap rmap,gmap,bmap,
       will generate a new 8-bit-per, dwide * dhigh image, which is dithered
       using colors found in rdmap, gdmap, bdmap arrays */

    /* returns ptr to a dwide*dhigh array of bytes, or NULL on failure */
    byte * picSmooth = Smooth(picSrc8, swide, shigh, dwide, dhigh);
    if (picSmooth)
    {
    	byte * picDithered = NULL;
    	DoColorDither(picSmooth, picDithered, dwide, dhigh);
        free(picSmooth);
        return picDithered;
    }

    return (byte *) NULL;
}



/***************************************************/
byte *Smooth(byte *picSrc8, int swide, int shigh, int dwide, int dhigh)
{
    /* does a SMOOTH resize from pic824 (which is either a swide*shigh, 8-bit
       pic, with colormap rmap,gmap,bmap OR a swide*shigh, 24-bit image, based
       on whether 'is24' is set) into a dwide * dhigh 24-bit image

       returns a dwide*dhigh 24bit image, or NULL on failure (malloc) */
    /* rmap,gmap,bmap should be 'desired' colors */

    byte *picSmooth, *ptrPicSmooth;
    int  *cxtab, *pxtab;
    size_t   y1Off, cyOff;
    size_t   ex, ey, cx, cy, px, py, apx, apy, x1, y1;
    size_t   cA, cB, cC, cD;
    size_t   pA, pB, pC, pD;
    int   retval;

    cA = cB = cC = cD = 0;
    size_t picSmoothSize=((size_t)dwide) * dhigh * 3;
    ptrPicSmooth = picSmooth = (byte *) malloc(picSmoothSize);
    if (!picSmooth)
    	REPORT_ERROR(ERR_MEM_NOTENOUGH,formatString("Unable to alloc: %lu",picSmoothSize));

    /* decide which smoothing routine to use based on type of expansion */
    if (dwide <  swide && dhigh <  shigh)
        retval = SmoothXY(picSmooth, picSrc8, swide, shigh, dwide, dhigh);

    else if (dwide <  swide && dhigh >= shigh)
        retval = SmoothX(picSmooth, picSrc8, swide, shigh, dwide, dhigh);

    else if (dwide >= swide && dhigh <  shigh)
        retval = SmoothY(picSmooth, picSrc8, swide, shigh, dwide, dhigh);

    else
    {
        /* dwide >= swide && dhigh >= shigh */

        /* cx,cy = original pixel in pic824.  px,py = relative position
           of pixel ex,ey inside of cx,cy as percentages +-50%, +-50%.
           0,0 = middle of pixel */

        /* we can save a lot of time by precomputing cxtab[] and pxtab[], both
           dwide arrays of ints that contain values for the equations:
             cx = (ex * swide) / dwide;
             px = ((ex * swide * 100) / dwide) - (cx * 100) - 50; */

        cxtab = (int *) malloc(dwide * sizeof(int));
        if (!cxtab)
        	REPORT_ERROR(ERR_MEM_NOTENOUGH,"Unable to alloc for smoothing");

        pxtab = (int *) malloc(dwide * sizeof(int));
        if (!pxtab)
        	REPORT_ERROR(ERR_MEM_NOTENOUGH,"Unable to alloc for smoothing");

        for (ex = 0; ex < dwide; ex++)
        {
            cxtab[ex] = (ex * swide) / dwide;
            pxtab[ex] = (((ex * swide) * 100) / dwide)
                        - (cxtab[ex] * 100) - 50;
        }

        for (ey = 0; ey < dhigh; ey++)
        {
            cy = (ey * shigh) / dhigh;
            py = (((ey * shigh) * 100) / dhigh) - (cy * 100) - 50;
            if (py < 0)
            {
                y1 = cy - 1;
                if (y1 < 0)
                    y1 = 0;
            }
            else
            {
                y1 = cy + 1;
                if (y1 > shigh - 1)
                    y1 = shigh - 1;
            }

            cyOff = (size_t) cy * swide ;    /* current line */
            y1Off = (size_t) y1 * swide ;    /* up or down one line, depending */

            /*      if ((ey&15) == 0) WaitCursor(); */

            for (ex = 0; ex < dwide; ex++)
            {
                byte *pptr, rA, gA, bA, rB, gB, bB, rC, gC, bC, rD, gD, bD;

                cx = cxtab[ex];
                px = pxtab[ex];

                if (px < 0)
                {
                    x1 = cx - 1;
                    if (x1 < 0)
                        x1 = 0;
                }
                else
                {
                    x1 = cx + 1;
                    if (x1 > swide - 1)
                        x1 = swide - 1;
                }

                cA = picSrc8[y1Off + x1];   /* corner pixel */
                cB = picSrc8[y1Off + cx];   /* up/down center pixel */
                cC = picSrc8[cyOff + x1];   /* left/right center pixel */
                cD = picSrc8[cyOff + cx];   /* center pixel */

                /* quick check */
                if (cA == cB && cB == cC && cC == cD)
                {
                    /* set this pixel to the same color as in pic8 */
                    *ptrPicSmooth++ = cD;
                    *ptrPicSmooth++ = cD;
                    *ptrPicSmooth++ = cD;
                }

                else
                {
                    /* compute weighting factors */
                    apx = abs(px);
                    apy = abs(py);
                    pA = (apx * apy) / 100;
                    pB = (apy * (100 - apx)) / 100;
                    pC = (apx * (100 - apy)) / 100;
                    pD = 100 - (pA + pB + pC);

                    byte val=(byte)((pA * cA) / 100 + (pB * cB) / 100 +
                                    (pC * cC) / 100 + (pD * cD) / 100);
                    /*
                    *ptrPicSmooth++ = (pA * rmap[cA]) / 100 + (pB * rmap[cB]) / 100 +
                    (pC * rmap[cC]) / 100 + (pD * rmap[cD]) / 100;

                    *pp++ = (pA * gmap[cA]) / 100 + (pB * gmap[cB]) / 100 +
                    (pC * gmap[cC]) / 100 + (pD * gmap[cD]) / 100;

                    *pp++ = (pA * bmap[cA]) / 100 + (pB * bmap[cB]) / 100 +
                    (pC * bmap[cC]) / 100 + (pD * bmap[cD]) / 100;                    }
                    */
                    *ptrPicSmooth++=val;
                    *ptrPicSmooth++=val;
                    *ptrPicSmooth++=val;
                }
            }
        }

        free(cxtab);
        free(pxtab);
        retval = 0;    /* okay */
    }

    if (retval)
    {    /* one of the Smooth**() methods failed */
        free(picSmooth);
        picSmooth = (byte *) NULL;
    }

    return picSmooth;
}

/***************************************************/
int SmoothX(byte *picSmooth, byte *picSrc8,
            int swide, int shigh, int dwide, int dhigh)
{
    byte *cptr, *cptr1;
    int  i, j;
    int  *lbufR, *lbufG, *lbufB;
    int  pixR, pixG, pixB;
    int  pcnt0, pcnt1, lastpix, pixcnt, thisline, ypcnt;
    int  *pixarr, *paptr;

    /* returns '0' if okay, '1' if failed (malloc) */

    /* for case where pic8 is shrunk horizontally and stretched vertically
       maps pic8 into an dwide * dhigh 24-bit picture.  Only works correctly
       when swide>=dwide and shigh<=dhigh */


    /* malloc some arrays */
    lbufR = (int *) calloc(swide, sizeof(int));
    lbufG = (int *) calloc(swide, sizeof(int));
    lbufB = (int *) calloc(swide, sizeof(int));
    pixarr = (int *) calloc(swide + 1, sizeof(int));
    if (!lbufR || !lbufG || !lbufB || !pixarr)
    	REPORT_ERROR(ERR_MEM_NOTENOUGH,"Cannot allocate memory for smoothing");

    for (j = 0; j <= swide; j++)
        pixarr[j] = (j * dwide + (15 * swide) / 16) / swide;

    cptr = picSrc8;
    cptr1 = cptr + swide;

    for (i = 0; i < dhigh; i++)
    {
        ypcnt = (((i * shigh) << 6) / dhigh) - 32;
        if (ypcnt < 0)
            ypcnt = 0;

        pcnt1 = ypcnt & 0x3f;                     /* 64ths of NEXT line to use */
        pcnt0 = 64 - pcnt1;                       /* 64ths of THIS line to use */

        thisline = ypcnt >> 6;

        cptr  = picSrc8 + thisline * swide;
        if (thisline + 1 < shigh)
            cptr1 = cptr + swide;
        else
            cptr1 = cptr;

        for (j = 0; j < swide; j++, cptr++, cptr1++)
        {
            lbufR[j] = ((*cptr * pcnt0) + (*cptr1 * pcnt1)) >> 6;
            lbufG[j] = ((*cptr * pcnt0) + (*cptr1 * pcnt1)) >> 6;
            lbufB[j] = ((*cptr * pcnt0) + (*cptr1 * pcnt1)) >> 6;
        }

        pixR = pixG = pixB = pixcnt = lastpix = 0;

        for (j = 0, paptr = pixarr; j <= swide; j++, paptr++)
        {
            if (*paptr != lastpix)
            {   /* write a pixel to pic24 */
                *picSmooth++ = pixR / pixcnt;
                *picSmooth++ = pixG / pixcnt;
                *picSmooth++ = pixB / pixcnt;
                lastpix = *paptr;
                pixR = pixG = pixB = pixcnt = 0;
            }

            if (j < swide)
            {
                pixR += lbufR[j];
                pixG += lbufG[j];
                pixB += lbufB[j];
                pixcnt++;
            }
        }
    }

    free(lbufR);
    free(lbufG);
    free(lbufB);
    free(pixarr);
    return 0;
}

/***************************************************/
int SmoothY(byte *picSmooth, byte *picSrc8, int swide,
            int shigh, int dwide, int dhigh)
{
    byte *clptr, *cptr, *cptr1;
    int  i, j;
    int  *lbufR, *lbufG, *lbufB, *pct0, *pct1, *cxarr, *cxptr;
    int  lastline, thisline, linecnt;
    int  retval;


    /* returns '0' if okay, '1' if failed (malloc) */

    /* for case where pic8 is shrunk vertically and stretched horizontally
       maps pic8 into a dwide * dhigh 24-bit picture.  Only works correctly
       when swide<=dwide and shigh>=dhigh */

    retval = 0;   /* no probs, yet... */

    lbufR = lbufG = lbufB = pct0 = pct1 = cxarr = NULL;
    lbufR = (int *) calloc(dwide, sizeof(int));
    lbufG = (int *) calloc(dwide, sizeof(int));
    lbufB = (int *) calloc(dwide, sizeof(int));
    pct0  = (int *) calloc(dwide, sizeof(int));
    pct1  = (int *) calloc(dwide, sizeof(int));
    cxarr = (int *) calloc(dwide, sizeof(int));

    if (!lbufR || !lbufG || !lbufB || !pct0 || ! pct1 || !cxarr)
    	REPORT_ERROR(ERR_MEM_NOTENOUGH,"Cannot allocate memory for smoothing");

    for (i = 0; i < dwide; i++)
    {                /* precompute some handy tables */
        int cx64;
        cx64 = (((i * swide) << 6) / dwide) - 32;
        if (cx64 < 0)
            cx64 = 0;
        pct1[i] = cx64 & 0x3f;
        pct0[i] = 64 - pct1[i];
        cxarr[i] = cx64 >> 6;
    }

    lastline = linecnt = 0;
    for (i = 0, clptr = picSrc8; i <= shigh; i++, clptr += swide)
    {
        /*    if ((i&15) == 0) WaitCursor();*/

        thisline = (i * dhigh + (15 * shigh) / 16) / shigh;

        if (thisline != lastline)
        {  /* copy a line to pic24 */
            for (j = 0; j < dwide; j++)
            {
                *picSmooth++ = lbufR[j] / linecnt;
                *picSmooth++ = lbufG[j] / linecnt;
                *picSmooth++ = lbufB[j] / linecnt;
            }

            xvbzero((char *) lbufR, dwide * sizeof(int));  /* clear out line bufs */
            xvbzero((char *) lbufG, dwide * sizeof(int));
            xvbzero((char *) lbufB, dwide * sizeof(int));
            linecnt = 0;
            lastline = thisline;
        }


        for (j = 0, cxptr = cxarr; j < dwide; j++, cxptr++)
        {
            cptr  = clptr + *cxptr;
            if (*cxptr < swide - 1)
                cptr1 = cptr + 1;
            else
                cptr1 = cptr;

            lbufR[j] += (((*cptr * pct0[j]) + (*cptr1 * pct1[j])) >> 6);
            lbufG[j] += (((*cptr * pct0[j]) + (*cptr1 * pct1[j])) >> 6);
            lbufB[j] += (((*cptr * pct0[j]) + (*cptr1 * pct1[j])) >> 6);
        }

        linecnt++;
    }

    return retval;
}

/***************************************************/
int SmoothXY(byte *pic24, byte *pic824,
             int swide, int shigh, int dwide, int dhigh)
{
    byte *cptr;
    int  i, j;
    int  *lbufR, *lbufG, *lbufB;
    int  pixR, pixG, pixB;
    int  lastline, thisline, lastpix, linecnt, pixcnt;
    int  *pixarr, *paptr;

    /* returns '0' if okay, '1' if failed (malloc) */

    /* shrinks pic8 into a dwide * dhigh 24-bit picture.  Only works correctly
       when swide>=dwide and shigh>=dhigh (ie, the picture is shrunk on both
       axes) */

    /* malloc some arrays */
    lbufR = (int *) calloc(swide, sizeof(int));
    lbufG = (int *) calloc(swide, sizeof(int));
    lbufB = (int *) calloc(swide, sizeof(int));
    pixarr = (int *) calloc(swide + 1, sizeof(int));
    if (!lbufR || !lbufG || !lbufB || !pixarr)
    	REPORT_ERROR(ERR_MEM_NOTENOUGH,"Cannot allocate memory for smoothing");

    for (j = 0; j <= swide; j++)
        pixarr[j] = (j * dwide + (15 * swide) / 16) / swide;

    lastline = linecnt = pixR = pixG = pixB = 0;
    cptr = pic824;

    for (i = 0; i <= shigh; i++)
    {
        /*    if ((i&15) == 0) WaitCursor();*/

        thisline = (i * dhigh + (15 * shigh) / 16) / shigh;

        if ((thisline != lastline))
        {      /* copy a line to pic24 */
            pixR = pixG = pixB = pixcnt = lastpix = 0;

            for (j = 0, paptr = pixarr; j <= swide; j++, paptr++)
            {
                if (*paptr != lastpix)
                {                 /* write a pixel to pic24 */
                    *pic24++ = (pixR / linecnt) / pixcnt;
                    *pic24++ = (pixG / linecnt) / pixcnt;
                    *pic24++ = (pixB / linecnt) / pixcnt;
                    lastpix = *paptr;
                    pixR = pixG = pixB = pixcnt = 0;
                }

                if (j < swide)
                {
                    pixR += lbufR[j];
                    pixG += lbufG[j];
                    pixB += lbufB[j];
                    pixcnt++;
                }
            }

            lastline = thisline;
            xvbzero((char *) lbufR, swide * sizeof(int));  /* clear out line bufs */
            xvbzero((char *) lbufG, swide * sizeof(int));
            xvbzero((char *) lbufB, swide * sizeof(int));
            linecnt = 0;
        }

        if (i < shigh)
        {
            for (j = 0; j < swide; j++, cptr++)
            {
                lbufR[j] += *cptr;
                lbufG[j] += *cptr;
                lbufB[j] += *cptr;
            }

            linecnt++;
        }
    }

    free(lbufR);
    free(lbufG);
    free(lbufB);
    free(pixarr);
    return 0;
}

/********************************************/
void DoColorDither(byte *picSmooth, byte *&picDithered, int w, int h)
{
    /* takes a 24 bit picture, of size w*h, dithers with the colors in
       rdisp, gdisp, bdisp (which have already been allocated),
       and generates an 8-bit w*h image, which it returns.
       ignores input value 'pic8'
       returns NULL on error

       note: the rdisp,gdisp,bdisp arrays should be the 'displayed' colors,
       not the 'desired' colors

       if picSrc is NULL, uses the passed-in pic8 (an 8-bit image) as
       the source, and the rmap,gmap,bmap arrays as the desired colors */

    byte *np, *ep;
    short *cache;
    int r2, g2, b2;
    int *thisline, *nextline, *thisptr, *nextptr, *tmpptr;
    int  i, j, rerr, gerr, berr, pwide3;
    int  imax, jmax;
    int key;
    long cnt1, cnt2;

    cnt1 = cnt2 = 0;
    pwide3 = w * 3;
    imax = h - 1;
    jmax = w - 1;
    ep = picSmooth;

    /* attempt to malloc things */
    picDithered = (byte *) malloc((size_t) w * h);
    cache = (short *) calloc(2 << 14, sizeof(short));
    thisline = (int *) malloc(pwide3 * sizeof(int));
    nextline = (int *) malloc(pwide3 * sizeof(int));
    if (!cache || !picDithered || !thisline || !nextline)
    	REPORT_ERROR(ERR_MEM_NOTENOUGH,"Cannot allocate memory for smoothing");

    np = picDithered;

    /* get first line of picture */

    if (picSmooth)
    {
        for (j = pwide3, tmpptr = nextline; j; j--, ep++)
            *tmpptr++ = (int) * ep;
    }
    else
    {
        for (j = w, tmpptr = nextline; j; j--, ep++)
        {
            *tmpptr++ = (int) *ep;
            *tmpptr++ = (int) *ep;
            *tmpptr++ = (int) *ep;
        }
    }


    for (i = 0; i < h; i++)
    {
        np = picDithered + i * w;
        /*    if ((i&15) == 0) WaitCursor();*/

        tmpptr = thisline;
        thisline = nextline;
        nextline = tmpptr;   /* swap */

        if (i != imax)
        {  /* get next line */
            if (!picSmooth)
                for (j = w, tmpptr = nextline; j; j--, ep++)
                {
                    *tmpptr++ = (int) *ep;
                    *tmpptr++ = (int) *ep;
                    *tmpptr++ = (int) *ep;
                }
            else
                for (j = pwide3, tmpptr = nextline; j; j--, ep++)
                    *tmpptr++ = (int) * ep;
        }

        /* dither a line, doing odd-lines right-to-left (serpentine) */
        thisptr = (i & 1) ? thisline + w * 3 - 3 : thisline;
        nextptr = (i & 1) ? nextline + w * 3 - 3 : nextline;
        if (i&1)
            np += w - 1;


        for (j = 0; j < w; j++)
        {
            int k, d, mind, closest;

            r2 = *thisptr++;
            g2 = *thisptr++;
            b2 = *thisptr++;
            if (i&1)
                thisptr -= 6;  /* move left */

            /* map r2,g2,b2 components (could be outside 0..255 range)
            into 0..255 range */

            if (r2 < 0 || g2 < 0 || b2 < 0)
            {   /* are there any negatives in RGB? */
                if (r2 < g2)
                {
                    if (r2 < b2)
                        k = 0;
                    else
                        k = 2;
                }
                else
                {
                    if (g2 < b2)
                        k = 1;
                    else
                        k = 2;
                }

                switch (k)
                {
                case 0:
                    g2 -= r2;
                    b2 -= r2;
                    d = (abs(r2) * 3) / 2;    /* RED */
                    r2 = 0;
                    g2 = (g2 > d) ? g2 - d : 0;
                    b2 = (b2 > d) ? b2 - d : 0;
                    break;

                case 1:
                    r2 -= g2;
                    b2 -= g2;
                    d = (abs(g2) * 3) / 2;    /* GREEN */
                    r2 = (r2 > d) ? r2 - d : 0;
                    g2 = 0;
                    b2 = (b2 > d) ? b2 - d : 0;
                    break;

                case 2:
                    r2 -= b2;
                    g2 -= b2;
                    d = (abs(b2) * 3) / 2;    /* BLUE */
                    r2 = (r2 > d) ? r2 - d : 0;
                    g2 = (g2 > d) ? g2 - d : 0;
                    b2 = 0;
                    break;
                }
            }

            if (r2 > 255 || g2 > 255 || b2 > 255)
            {   /* any overflows in RGB? */
                if (r2 > g2)
                {
                    if (r2 > b2)
                        k = 0;
                    else
                        k = 2;
                }
                else
                {
                    if (g2 > b2)
                        k = 1;
                    else
                        k = 2;
                }

                switch (k)
                {
                case 0:
                    g2 = (g2 * 255) / r2;
                    b2 = (b2 * 255) / r2;
                    r2 = 255;
                    break;
                case 1:
                    r2 = (r2 * 255) / g2;
                    b2 = (b2 * 255) / g2;
                    g2 = 255;
                    break;
                case 2:
                    r2 = (r2 * 255) / b2;
                    g2 = (g2 * 255) / b2;
                    b2 = 255;
                    break;
                }
            }

            key = ((r2 & 0xf8) << 6) | ((g2 & 0xf8) << 1) | (b2 >> 4);
            if (key >= (2 << 14))
            {
                fprintf(stderr, "'key' overflow in DoColorDither()");
                exit(-1);
            }

            if (cache[key])
            {
                *np = (byte)(cache[key] - 1);
                cnt1++;
            }
            else
            {
                /* not in cache, have to search the colortable */
                cnt2++;

                mind = 10000;
                for (k = closest = 0; k < 256 && mind > 7; k++)
                {
                    d = abs(r2 - k)
                        + abs(g2 - k)
                        + abs(b2 - k);
                    if (d < mind)
                    {
                        mind = d;
                        closest = k;
                    }
                }
                cache[key] = closest + 1;
                *np = closest;
            }


            /* propogate the error */
            rerr = r2 - *np;
            gerr = g2 - *np;
            berr = b2 - *np;

            if (j != jmax)
            {  /* adjust LEFT/RIGHT pixel */
                thisptr[0] += (rerr / 2);
                rerr -= (rerr / 2);
                thisptr[1] += (gerr / 2);
                gerr -= (gerr / 2);
                thisptr[2] += (berr / 2);
                berr -= (berr / 2);
            }

            if (i != imax)
            { /* adjust BOTTOM pixel */
                nextptr[0] += rerr;    /* possibly all err if we're at l/r edge */
                nextptr[1] += gerr;
                nextptr[2] += berr;
            }

            if (i&1)
            {
                nextptr -= 3;
                np--;
            }
            else
            {
                nextptr += 3;
                np++;
            }
        }
    }


    free(thisline);
    free(nextline);
    free(cache);
}
