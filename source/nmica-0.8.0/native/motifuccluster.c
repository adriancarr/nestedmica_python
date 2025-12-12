/*
 *                      NestedMICA
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <alloca.h>

#include "nmica.h"
#include "fastlog.h"

#include "net_derkholm_nmica_model_motif_MotifUncountedLikelihoodNative.h"

#define MICANATIVE_VERSION 34

#define M_SET(m, cols, row, col, val) m[((row) * (cols)) + (col)]=val
#define M_GET(m, cols, row, col) (m[((row) * (cols)) + (col)])

#define LOG_2 0.69314718056
#define ONE_LOG_2 1.44269504089
#define NLOG2 fastlog2

static inline jdouble addLog2(jdouble x, jdouble y) {    
        jdouble base;
    
        if (x <= NMICA_NEGATIVE_INFINITY) {
            return y;
        } else if (y <= NMICA_NEGATIVE_INFINITY) {
            return x;
        } else {
            if (x > y) {
	        base = x;
            } else {
                base = y;
            }
            return NLOG2(exp((x - base) * LOG_2) + exp((y - base) * LOG_2)) + base;
        }
}

static inline jdouble multiAddLog2(jdouble *terms, int n) {    
        jdouble base = NMICA_NEGATIVE_INFINITY;
        jdouble tot = 0;
        int i;
        
        for (i = 0; i < n; ++i) {
            if (terms[i] > base) {
                base = terms[i];
            }
        }
        
        for (i = 0; i < n; ++i) {
            if (terms[i] > NMICA_NEGATIVE_INFINITY) {
                tot += exp((terms[i] - base) * LOG_2);
            }
        }
        
        return NLOG2(tot) + base;
}

static void sumMotifsMulti(
	  jint length,
      jdouble *bgScores, 
      jint numMotifs,
      jdouble *wmScores, 
      jint wmScores_rows, 
      jint wmScores_columns, 
      jint *advances, 
      jdouble *matrix,
      jint matrix_cols,
      jdouble clusterInPenalty,
      jdouble clusterOutPenalty,
      jdouble outsideBasePenalty,
      jdouble insideBasePenalty,
      jdouble *motifPenalty
  )
{
    jint i, m, wml;
    jdouble fromScore, emitScore;
    jdouble *scoreTerms = alloca((numMotifs + 2) * sizeof(jdouble));
    jdouble bgEmit;
    
    for (i = 1; i <= length; ++i) {
    	bgEmit = bgScores[i - 1];
    	
    	/* inside */
        scoreTerms[numMotifs + 1] = M_GET(matrix, matrix_cols, 0, i - 1) + bgEmit+ insideBasePenalty;
        scoreTerms[numMotifs] = M_GET(matrix, matrix_cols, 1, i - 1) + bgEmit + clusterInPenalty;
        for (m = 0; m < numMotifs; ++m) {
            wml = advances[m];
            if (i >= wml) {
                fromScore = matrix[i - wml];
                emitScore = M_GET(wmScores, wmScores_columns, i - 1, m);
                scoreTerms[m] = fromScore + emitScore + motifPenalty[m];
            } else {
                scoreTerms[m] = NMICA_NEGATIVE_INFINITY;
            }
        }
        M_SET(matrix, matrix_cols, 0, i, multiAddLog2(scoreTerms, numMotifs + 2));
        
        /* outside */
        
        scoreTerms[0] = M_GET(matrix, matrix_cols, 1, i - 1) + bgEmit + outsideBasePenalty;
        scoreTerms[1] = M_GET(matrix, matrix_cols, 0, i - 1) + bgEmit + clusterOutPenalty;
        M_SET(matrix, matrix_cols, 1, i, multiAddLog2(scoreTerms, 2));
    }
}

JNIEXPORT jdouble JNICALL 
  Java_net_derkholm_nmica_model_motif_MotifUncountedClusterLikelihoodNative_nativeSumMotifs(
      JNIEnv *env,
      jobject thiz,
      jint length,
      jdoubleArray bgScoresJ, 
      jint numMotifs,
      jdoubleArray wmScoresJ, 
      jint wmScores_rows, 
      jint wmScores_columns, 
      jintArray advancesJ, 
      jdoubleArray motifTransJ,
      jdouble clusterIn,
      jdouble clusterOut
  )
{
    jint *advances;
    jdouble *wmScores;
    jdouble *bgScores;
    jdouble *matrix;
    jdouble *motifTrans;
    
    jint m;
    jdouble sumTrans = 0;
    jdouble* motifPenalty;
    jdouble clusterInPenalty, clusterOutPenalty;
    jdouble insideBasePenalty;
    jdouble outsideBasePenalty;
    jdouble hood;
    
    /* Pin the arrays we need */
    
    advances = (*env)->GetIntArrayElements(
        env,
        advancesJ,
        NULL
    );
    wmScores = (*env)->GetDoubleArrayElements(
        env,
        wmScoresJ,
        NULL
    );
    bgScores = (*env)->GetDoubleArrayElements(
        env,
        bgScoresJ,
        NULL
    );
    motifTrans = (*env)->GetDoubleArrayElements(
        env,
        motifTransJ,
        NULL
    );
    
    motifPenalty = alloca(numMotifs * sizeof(jdouble));
    for (m = 0; m < numMotifs; ++m) {
        motifPenalty[m] = log(motifTrans[m]) * ONE_LOG_2;
        sumTrans += motifTrans[m];
    }
    clusterInPenalty = log(clusterIn) * ONE_LOG_2;
    clusterOutPenalty = log(clusterOut) * ONE_LOG_2;
    outsideBasePenalty = log(1.0 - clusterIn) * ONE_LOG_2;
    insideBasePenalty = log(1.0 - sumTrans - clusterOut) * ONE_LOG_2;
    
    int matrix_cols = length + 1;
    matrix = alloca(matrix_cols * sizeof(jdouble) * 2);
    M_SET(matrix, matrix_cols, 0, 0, log(0.5) * ONE_LOG_2);
    M_SET(matrix, matrix_cols, 1, 0, log(0.5) * ONE_LOG_2);
    
    sumMotifsMulti(length, bgScores, numMotifs, wmScores, wmScores_rows, wmScores_columns, advances, matrix, matrix_cols, clusterInPenalty, clusterOutPenalty, outsideBasePenalty, insideBasePenalty, motifPenalty);
    
    hood = addLog2(M_GET(matrix, matrix_cols, 0, length), M_GET(matrix, matrix_cols, 1, length));
    
    /* Undo the pinning */
    (*env)->ReleaseDoubleArrayElements(env, bgScoresJ, bgScores, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, wmScoresJ, wmScores, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, advancesJ, advances, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, motifTransJ, motifTrans, JNI_ABORT);
    
    return hood;
}

