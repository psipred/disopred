/* Neural Network Identification of Protein Disorder */
/* Copyright (C) 1998 David T. Jones - Last Edited : Nov 2014 */

/* Combined prediction module */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "disordcomb_net.h"

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

enum
{
    FALSE, TRUE
};

#define MAXSQLEN 50000

#define SQR(x) ((x)*(x))

#define REAL float

void           *calloc(), *malloc();

char           *wtfnm;

int             nwtsum, fwt_to[TOTAL], lwt_to[TOTAL];
REAL            blrate = ILRATE, lrate = ILRATE, alpha = IALPHA;
REAL            activation[TOTAL], bias[TOTAL], netinput[TOTAL], *weight[TOTAL];
REAL            bias[TOTAL];

int             nhelix, nsheet, ncoil, restot, nblocks;

float           pred1[MAXSQLEN], pred2[MAXSQLEN], pred3[MAXSQLEN];

char            seq[MAXSQLEN], ssstruc[MAXSQLEN], ssstrel[MAXSQLEN];
int             seqlen;

char           *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "UNK"
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};

void fail(char *s)
{
    fprintf(stderr, "%s\n", s);
    exit(1);
}


/*
 * Back Propagation Engine - as described by McClelland & Rumelhart /
 * Sejnowski & Rosenberg
 */

/* logistic 'squashing' function (+/- 1.0) */
#define logistic(x) ((REAL)1.0 / ((REAL)1.0 + (REAL)exp(-(x))))

void
                compute_output()
{
    int             i, j;
    REAL            netinp;

    for (i = NUM_IN; i < TOTAL; i++)
    {
	netinp = bias[i];

	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	    netinp += activation[j] * weight[i][j];

	/* Trigger neuron */
	activation[i] = logistic(netinp);
    }
}

/*
 * load weights - load all link weights from a disk file
 */
void                load_wts(char *fname)
{
    int             i, j;
    double          t;
    FILE           *ifp;

    if (!(ifp = fopen(fname, "r")))
	fail("Cannot open weights file!");

    /* Load input units to hidden layer weights */
    for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	{
	    if (fscanf(ifp, "%lf", &t) != 1)
		fail("Bad weights file!");
	    weight[i][j] = t;
	}

    /* Load hidden layer to output units weights */
    for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	{
	    if (fscanf(ifp, "%lf", &t) != 1)
		fail("Bad weights file!");
	    weight[i][j] = t;
	}

    /* Load bias weights */
    for (j = NUM_IN; j < TOTAL; j++)
    {
	if (fscanf(ifp, "%lf", &t) != 1)
	    fail("Bad weights file!");
	bias[j] = t;
    }

    fclose(ifp);
}

/* Initialize network - wire up units and make random weights */
void                init()
{
    int             i;

    for (i = NUM_IN; i < TOTAL; i++)
	if (!(weight[i] = calloc(TOTAL - NUM_OUT, sizeof(REAL))))
	    fail("init: Out of Memory!");

    /* Connect input units to hidden layer */
    for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
    {
	fwt_to[i] = 0;
	lwt_to[i] = NUM_IN;
    }

    /* Connect hidden units to output layer */
    for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
    {
	fwt_to[i] = NUM_IN;
	lwt_to[i] = NUM_IN + NUM_HID;
    }
}

/* Convert AA letter to numeric code (0-20) */
int                aanum(int ch)
{
    static int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
	20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}


/* Make prediction */
void                predict()
{
    int             j, winpos;
    char            pred;
    float           prec;

    printf("#         ----- DISOPRED version 3.1 -----\n");
    printf("# Disordered residues are marked with asterisks (*)\n");
    printf("#    Ordered residues are marked with dots (.)\n");

    for (winpos = 0; winpos < seqlen; winpos++)
    {
	for (j = 0; j < NUM_IN; j++)
	    activation[j] = 0.0;
	for (j = WINL; j <= WINR; j++)
	{
	    if (j + winpos >= 0 && j + winpos < seqlen)
	    {
		activation[(j - WINL) * IPERGRP] = pred1[j + winpos];
		activation[(j - WINL) * IPERGRP+1] = pred2[j + winpos];
		activation[(j - WINL) * IPERGRP+2] = pred3[j + winpos];
	    }
	    else
		activation[(j - WINL) * IPERGRP + 3] = 1.0;
	}
	compute_output();

	prec = 0.5 + 0.5 * (activation[TOTAL - NUM_OUT + 1] - activation[TOTAL - NUM_OUT]);
	
	if (prec < 0.5)
	    pred = '.';
	else
	    pred = '*';

	printf("%5d %c %c %4.2f\n", winpos+1, seq[winpos], pred, prec);
    }
}

int main(int argc, char **argv)
{
    int  i;
    char buf[512];
    FILE *ifp;

    /* malloc_debug(3); */
    if (argc != 5)
	fail("usage : disordcomb_pred weight-file disopred2file nndisopredfile neighbfile");
    init();
    load_wts(wtfnm = argv[1]);

    if (!(ifp = fopen(argv[2], "r")))
	fail("Cannot open diso file!");
    i = 0;
    while (!feof(ifp))
    {
	if (!fgets(buf, 512, ifp))
	    break;
	
	if (buf[0] != ' ' || !isdigit(buf[4]))
	    continue;
	
	if (sscanf(buf, "%*s%*s%*s%f", &pred1[i]) != 1)
	    fail("Bad diso file!");
	
	seq[i++] = buf[6];
    }
    fclose(ifp);
    
    seqlen = i;

    if (!(ifp = fopen(argv[3], "r")))
	fail("Cannot open nndiso file!");
    i = 0;
    while (!feof(ifp))
    {
	if (!fgets(buf, 512, ifp))
	    break;
	
	if (buf[0] != ' ' || !isdigit(buf[4]))
	    continue;
	
	if (sscanf(buf, "%*s%*s%*s%f", &pred2[i]) != 1)
	    fail("Bad nndiso file!");
	
	if (buf[8] == '.')
	    pred2[i] = -pred2[i];
	
	i++;
    }
    fclose(ifp);

    if (i != seqlen)
	fail("nndiso file length mismatch!");

    if (!(ifp = fopen(argv[4], "r")))
	fail("Cannot open dnb file!");
    i = 0;
    while (!feof(ifp))
    {
	if (!fgets(buf, 512, ifp))
	    break;
	
	if (buf[0] != ' ' || !isdigit(buf[4]))
	    continue;

	if (sscanf(buf, "%*s%*s%*s%f", &pred3[i]) != 1)
	    fail("Bad dnb file!");
	
	pred3[i] *= 0.001;

	if (buf[8] == '.')
	    pred3[i] = -pred3[i];
	    
	i++;
    }
    fclose(ifp);

    if (i != seqlen)
	fail("dnb file length mismatch!");

    predict();

    return 0;
}
