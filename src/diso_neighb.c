/* Neural Network Identification of Secondary Structure */
/* Copyright (C) 1998 David T. Jones - Created : March 1998 */

/* Version 3.0 - Created : Feb 2010 */

/* Training Module */

/*
 * Description: This program provides an implementation of a neural network
 * containing one hidden layer, which uses the generalized backpropagation
 * delta rule for learning.
 */

/* This version is hardcoded as a two-layer NN to allow code vectorisation */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

enum
{
    FALSE, TRUE
};

#define MAXSEQLEN 50000
#define DBSIZE 20000

#define WLEN 7

#define BIG 999999

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SQR(x) ((x)*(x))

/* logistic 'squashing' function (+/- 1.0) */
#define logistic(x) (1.0 / (1.0 + exp(-(x))))

void           *calloc(), *malloc();

int nblocks, restot, seqlen;

int             profile[MAXSEQLEN][20];

char seq[MAXSEQLEN];


struct entry
{
    char           *id, *seq, *disord;
    int           **profile, length;
}
db[DBSIZE];

char           *rnames[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "???"
};

/*  BLOSUM 62 */
const short           aamat[23][23] =
{
    {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0},
    {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1},
    {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1},
    {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1},
    {0, -3, -3, -3,10, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2},
    {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1},
    {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1},
    {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1},
    {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1},
    {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1},
    {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1,  3, -1, -3, -3, -1},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2},
    {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0},
    {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2},
    {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1},
    {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1},
    {-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1},
    {-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1},
    {0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, 4}
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};

const char *rescodes = "ARNDCQEGHILKMFPSTWYVXXX";

void            err(char *s)
{
    fprintf(stderr, "%s\n", s);
}

void            fail(char *s)
{
    fprintf(stderr, "%s\n", s);
    exit(1);
}

/* Convert AA letter to numeric code (0-20) */
int
aanum(ch)
    int             ch;
{
    static const int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
	20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Make prediction */
void            predict()
{
    int             i, j, k, n;
    int sim, maxsim, best_n = 0, best_j = 0;
    int pscore_d[MAXSEQLEN], pscore_o[MAXSEQLEN];
    
    for (i=0; i<seqlen; i++)
	pscore_d[i] = pscore_o[i] = -BIG;
    
    for (i=0; i<=seqlen-WLEN; i++)
    {
	maxsim = -BIG;

	for (n = 0; n<nblocks; n++)
	{
	    for (j=0; j<=db[n].length-WLEN; j++)
	    {
		if (i < 10 && j > 10)
		    continue;
		if (i > seqlen-WLEN-10 && j < db[n].length-WLEN-10)
		    continue;
		
		sim = 0;
		for (k=0; k<WLEN; k++)
		{
		    if (seq[i+k] < 20)
			sim += db[n].profile[j+k][(int)seq[i+k]];
		    if (db[n].seq[j+k] < 20)
			sim += profile[i+k][(int)db[n].seq[j+k]];
		}
		sim /= 2*WLEN;

		if (sim > maxsim)
		{
		    maxsim = sim;
		    best_n = n;
		    best_j = j;
		}
	    }
	}

	j = best_j;
	n = best_n;
	sim = maxsim;
	
	for (k=0; k<WLEN; k++)
	{
	    if (db[n].disord[j+k] && sim > pscore_d[i+k])
		pscore_d[i+k] = sim;
	    if (!db[n].disord[j+k] && sim > pscore_o[i+k])
		pscore_o[i+k] = sim;
	}
    }

    for (i=0; i<seqlen; i++)
	printf("%5d %c %c %7d\n", i + 1, rescodes[(int)seq[i]], pscore_d[i] > pscore_o[i] ? '*' : '.', MAX(pscore_d[i], pscore_o[i]));
}

#define CH malloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);

void           *allocmat(int rows, int columns, int size, int clrflg)
{
    int             i;
    void          **p;

    p = malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    if (clrflg)
    {
	for (i = 0; i < rows; i++)
	    if ((p[i] = calloc(columns, size)) == NULL)
		fail("allocmat: calloc [][] failed!");
    }
    else
	for (i = 0; i < rows; i++)
	    if ((p[i] = malloc(columns * size)) == NULL)
		fail("allocmat: malloc [][] failed!");

    return p;
}

/* Read in data */
void            read_dat(char *filename)
{
    FILE           *ifp, *tfp;
    char            dsoid[20], buf[MAXSEQLEN], tdbname[512], *tdbpath;
    int             i, j, nres, consv[20];
    
    if (!(ifp = fopen(filename, "r")))
	fail("Cannot open training list file!");

    if (!(tdbpath = getenv("DSO_LIB_PATH")))
        fail("Cannot find the dso_lib path!");


    while (!feof(ifp))
    {
	if (fscanf(ifp, "%s", dsoid) != 1)
	    break;
	if (dsoid[0] == '#')
	    continue;
	strcpy(tdbname, tdbpath);
//	strcpy(tdbname, "./data/");
	strcat(tdbname, dsoid);
	strcat(tdbname, ".dso");

/*	puts(tdbname); */
	if (!(tfp = fopen(tdbname, "r")))
	{
	    puts(tdbname);
	    fail("Cannot open DSO file!");
	}
	if (!fgets(buf, 512, tfp))
	    break;
	sscanf(buf, "%*s%*s%d", &nres);
	if (!(db[nblocks].id = strdup(dsoid)))
	    fail("Out of memory!");
	db[nblocks].length = nres;
	if (!(db[nblocks].seq = malloc(nres))
	    || !(db[nblocks].disord = malloc(nres)))
	    fail("Out of memory!");
	db[nblocks].profile = allocmat(nres, 20, sizeof(int), FALSE);

	i = 0;
	while (!feof(tfp))
	{
	    if (!fgets(buf, 512, tfp))
		break;
	    if (sscanf(buf + 15, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
		       &consv[0], &consv[1], &consv[2], &consv[3], &consv[4],
		       &consv[5], &consv[6], &consv[7], &consv[8], &consv[9],
		       &consv[10], &consv[11], &consv[12], &consv[13],
		       &consv[14], &consv[15], &consv[16], &consv[17],
		       &consv[18], &consv[19]) != 20)
	    {
	        puts(dsoid);
		fail("Bad consensus records in dso file!");
	    }
	    else
	    {
		for (j = 0; j < 20; j++)
		    db[nblocks].profile[i][j] = consv[j];
	    }

	    if (buf[5] == '-')
	    {
		db[nblocks].seq[i] = 20;
		db[nblocks].disord[i] = 0;
	    }
	    else
	    {
		db[nblocks].seq[i] = j = aanum(buf[5]);
		db[nblocks].disord[i] = buf[7] != '.';
	    }
	    i++;
	}
	fclose(tfp);

	if (i != nres)
	    continue;

	restot += nres;
	nblocks++;
    }

    fclose(ifp);
}


/* Read PSI AA frequency data */
int             getmtx(FILE *lfil)
{
    int             j, naa;
    char            buf[256];

    if (fscanf(lfil, "%d", &naa) != 1)
      fail("Bad mtx file - no sequence length!");

    if (naa > MAXSEQLEN)
      fail("Input sequence too long!");

    if (fscanf(lfil, "%s", seq) != 1)
      fail("Bad mtx file - no sequence!");

    while (!feof(lfil))
      {
	if (!fgets(buf, 65536, lfil))
	  fail("Bad mtx file!");
	if (!strncmp(buf, "-32768 ", 7))
	  {
	    for (j=0; j<naa; j++)
	      {
		if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &profile[j][ALA],  &profile[j][CYS], &profile[j][ASP],  &profile[j][GLU],  &profile[j][PHE],  &profile[j][GLY],  &profile[j][HIS],  &profile[j][ILE],  &profile[j][LYS],  &profile[j][LEU],  &profile[j][MET],  &profile[j][ASN],  &profile[j][PRO],  &profile[j][GLN],  &profile[j][ARG],  &profile[j][SER],  &profile[j][THR],  &profile[j][VAL],  &profile[j][TRP],  &profile[j][TYR]) != 20)
		  fail("Bad mtx format!");
		seq[j] = aanum(seq[j]);
		if (!fgets(buf, 65536, lfil))
		  break;
	      }
	  }
      }

    return naa;
}


int main(int argc, char **argv)
{
    FILE *ifp;

    if (argc < 2)
	fail("usage : diso_neighb mtxfile  {dso.lst}");

    ifp = fopen(argv[1], "r");
    if (!ifp)
	exit(1);
    seqlen = getmtx(ifp);
    fclose(ifp);

    if (argc > 2)
	read_dat(argv[2]);
    else
	read_dat("dso.lst");

    predict();

    return 0;
}
