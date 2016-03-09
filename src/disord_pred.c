/* Neural Network Identification of Protein Disorder */
/* Copyright (C) 1998 David T. Jones - Last Edited : Nov 2014 */

/* Prediction Module */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "disord_net.h"

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

int             profile[MAXSQLEN][20];
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

void                err(char *s)
{
    fprintf(stderr, "%s\n", s);
}

void                fail(char *s)
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
void
                load_wts(char *fname)
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
void
                init()
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
int
                aanum(ch)
     int             ch;
{
    static int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
	20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}


/* Make prediction */
void
                predict()
{
    int             aa, j, winpos;
    char            pred;
    float           aacomp[20];

    for (aa = 0; aa < 20; aa++)
	aacomp[aa] = 0;
    for (winpos = 0; winpos < seqlen; winpos++)
	if (seq[winpos] < 20)
	    aacomp[(int)seq[winpos]]++;
    for (aa = 0; aa < 20; aa++)
	aacomp[aa] /= seqlen;
    
    for (winpos = 0; winpos < seqlen; winpos++)
    {
	for (j = 0; j < NUM_IN; j++)
	    activation[j] = 0.0;
	for (j = WINL; j <= WINR; j++)
	{
	    if (j + winpos >= 0 && j + winpos < seqlen)
	    {
		for (aa = 0; aa < 20; aa++)
		    activation[(j - WINL) * 21 + aa] = profile[j + winpos][aa] / 1000.0;
	    }
	    else
		activation[(j - WINL) * IPERGRP + 20] = 1.0;
	}
	for (aa = 0; aa < 20; aa++)
	    activation[(WINR - WINL + 1) * IPERGRP + aa] = aacomp[aa];
	activation[(WINR - WINL + 1) * IPERGRP + 20] = log(1.0 + winpos);
	activation[(WINR - WINL + 1) * IPERGRP + 21] = log(seqlen - winpos);
	compute_output();
	if (activation[TOTAL - NUM_OUT] > activation[TOTAL - NUM_OUT + 1])
	    pred = '.';
	else
	    pred = '*';

	printf("%5d %c %c %4.2f\n", winpos+1, seq[winpos], pred, fabs(activation[TOTAL - NUM_OUT] - activation[TOTAL - NUM_OUT + 1]));
    }
}

#define CH malloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);

/* Return current CPU time usage - extremely portable version */

#ifdef unix
#define CLKRATE 1000000
#endif

#ifdef CLK_TCK
#define CLKRATE (CLK_TCK)
#endif

double
                cpu_time()
{
    static double   runtime = 0.0;
    static unsigned long last_t;
    unsigned long   temp;

    temp = clock() & 0x7fffffff;
    if (temp >= last_t)
	runtime += (double) abs(temp - last_t) / CLKRATE;
    else
	runtime += (double) (0x7fffffff - last_t + (temp + 1)) / CLKRATE;
    last_t = temp;
    return (runtime);
}

/* Read NBRF formatted sequence data */
void
                getseq(char *seq, FILE * ifp)
{
    char            ch;

    while (!feof(ifp) && (ch = getc(ifp)) != '*')
	if (isupper(ch))
	    *seq++ = ch;
    (void) getc(ifp);
    *seq++ = '\0';
}

/* Read PSI AA frequency data */
int             getmtx(FILE *lfil)
{
    int             j, naa;
    char            buf[512];
    
    if (fscanf(lfil, "%d", &naa) != 1)
	fail("Bad mtx file!");
    
    if (naa > MAXSQLEN)
	fail("Input sequence too long!");
    
    if (fscanf(lfil, "%s", seq) != 1)
	fail("Bad mtx file!");
    
    while (!feof(lfil))
    {
	if (!fgets(buf, 512, lfil))
	    fail("Bad mtx file!");
	if (!strncmp(buf, "-32768 ", 7))
	{
	    for (j=0; j<naa; j++)
	    {
		if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &profile[j][ALA],  &profile[j][CYS], &profile[j][ASP],  &profile[j][GLU],  &profile[j][PHE],  &profile[j][GLY],  &profile[j][HIS],  &profile[j][ILE],  &profile[j][LYS],  &profile[j][LEU],  &profile[j][MET],  &profile[j][ASN],  &profile[j][PRO],  &profile[j][GLN],  &profile[j][ARG],  &profile[j][SER],  &profile[j][THR],  &profile[j][VAL],  &profile[j][TRP],  &profile[j][TYR]) != 20)
		    fail("Bad PSI BLAST output!");
		if (!fgets(buf, 512, lfil))
		    break;
	    }
	}
    }
    
    return naa;
}


/* Read PSIPRED VFORMAT prediction data */
int             getpsipredv(FILE * lfil)
{
    int             naa;
    float confc, confh, confe;
    char            buf[256];

    if (!fgets(buf, 256, lfil))
	fail("Bad PSIPRED VFORMAT file!");

    if (!fgets(buf, 256, lfil))
	fail("Bad PSIPRED VFORMAT file!");

    naa = 0;
    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    break;
	if (sscanf(buf+10, "%f%f%f", &confc, &confh, &confe) != 3)
	    break;
	switch (buf[7])
	{
	case 'H':
	    ssstruc[naa] = 'H';
	    break;
	case 'E':
	    ssstruc[naa] = 'E';
	    break;
	default:
	    ssstruc[naa] = 'C';
	    break;
	}
	ssstrel[naa++] = 100*(2*MAX(MAX(confc, confh),confe)-(confc+confh+confe)-MIN(MIN(confc, confh),confe));
    }

    return naa;
}

int main(int argc, char **argv)
{
    FILE           *ifp;

    /* malloc_debug(3); */
    if (argc != 3)
	fail("usage : diso_pred weight-file mtxfile");

    init();
    load_wts(wtfnm = argv[1]);

    ifp = fopen(argv[2], "r");
    if (!ifp)
	fail("Cannot open seq file!");
    seqlen = getmtx(ifp);

    fclose(ifp);
    
    predict();

    return 0;
}
