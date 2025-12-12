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

/*

jint JNI_OnLoad(JavaVM *vm, void *reserved) {
    fprintf(stderr, "Loading MotifUncountedLikelihood native DP version %d\n", MICANATIVE_VERSION);
    return JNI_VERSION_1_4;
}

*/

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

/*

static void sumMotifsSingle(
	  jint length,
      jdouble *bgScores, 
      jint numMotifs,
      jdouble *wmScores, 
      jint wmScores_rows, 
      jint wmScores_columns, 
      jint *advances, 
      jdouble *matrix,
      jdouble basePenalty,
      jdouble motifPenalty
  )
{
    int wml = advances[0];
    int i;
    jdouble score;
    jdouble fromScore, emitScore;
    
    for (i = 1; i <= length; ++i) {
        score = matrix[i - 1] + bgScores[i - 1] + basePenalty;
        if (i >= wml) {
            fromScore = matrix[i - wml];
            emitScore = wmScores[i - 1];
            score = addLog2(score, fromScore + emitScore + motifPenalty);
        }
        matrix[i] = score;
    }
}

*/

static void sumMotifsMulti(
	  jint length,
      jdouble *bgScores, 
      jint numMotifs,
      jdouble *wmScores, 
      jint wmScores_rows, 
      jint wmScores_columns, 
      jint *advances, 
      jdouble *matrix,
      jdouble basePenalty,
      jdouble *motifPenalty
  )
{
    jint i, m, wml;
    jdouble fromScore, emitScore;
    jdouble *scoreTerms = alloca((numMotifs + 1) * sizeof(jdouble));
    
    for (i = 1; i <= length; ++i) {
        scoreTerms[numMotifs] = matrix[i - 1] + bgScores[i - 1] + basePenalty;
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
        matrix[i] = multiAddLog2(scoreTerms, numMotifs + 1);
    }
}

JNIEXPORT jdouble JNICALL 
  Java_net_derkholm_nmica_model_motif_MotifUncountedLikelihoodNative_nativeSumMotifs(
      JNIEnv *env,
      jobject thiz,
      jint length,
      jdoubleArray bgScoresJ, 
      jint numMotifs,
      jdoubleArray wmScoresJ, 
      jint wmScores_rows, 
      jint wmScores_columns, 
      jintArray advancesJ, 
      jdoubleArray motifTransJ
  )
{
    jint stateSpace;
    jint* advances;
    jdouble* wmScores;
    jdouble* bgScores;
    jdouble* motifTrans;
    jdouble* matrix;
    
    jint m;
    jdouble* motifPenalty;
    jdouble sumTrans = 0;
    jdouble basePenalty;
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
    basePenalty = log(1.0 - sumTrans) * ONE_LOG_2;
    
    matrix = alloca((length + 1) * sizeof(jdouble));
    matrix[0] = 0.0;
    
    sumMotifsMulti(length, bgScores, numMotifs, wmScores, wmScores_rows, wmScores_columns, advances, matrix, basePenalty, motifPenalty);
    
    hood = matrix[length];
    
    /* Undo the pinning */
    (*env)->ReleaseDoubleArrayElements(env, bgScoresJ, bgScores, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, wmScoresJ, wmScores, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, advancesJ, advances, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, motifTransJ, motifTrans, JNI_ABORT);
    
    return hood;
}

