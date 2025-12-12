#include <stdio.h>
#include <math.h>

#include "nmica.h"
#include "fastlog.h"
#include "net_derkholm_nmica_maths_NativeMath.h"

#define LOG_2 0.69314718056
#define ONE_LOG_2 1.44269504089

jint JNI_OnLoad(JavaVM *vm, void *reserved) {
    /* fprintf(stderr, "Loading nmica native code\n"); */
    init_fastlog();
    return JNI_VERSION_1_2;
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_log
  (JNIEnv *env, jclass clazz, jdouble x)
{
    return log(x);
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_exp
  (JNIEnv *env, jclass clazz, jdouble x)
{
    return exp(x);
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_loopLog
  (JNIEnv *env, jclass clazz, jdouble x, jint iter)
{
    double dummy;
    int i;
    
    for (i = 0; i < iter; ++i) {
        dummy = log(x);
    }
    return dummy;
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_loopExp
  (JNIEnv *env, jclass clazz, jdouble x, jint iter)
{
    double dummy;
    int i;
    
    for (i = 0; i < iter; ++i) {
        dummy = exp(x);
    }
    return dummy;
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_addLog
  (JNIEnv *env, jclass clazz, jdouble x, jdouble y)
{
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
        return log(exp(x - base) + exp(y - base)) + base;
    }
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_log2
  (JNIEnv *env, jclass clazz, jdouble x)
{
  return log(x) * ONE_LOG_2;
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_exp2
  (JNIEnv *env, jclass clazz, jdouble x)
{
  return exp(x * LOG_2);
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_fastlog2
  (JNIEnv *env, jclass clazz, jdouble x)
{
  return fastlog2(x);
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_addLog2__DD
  (JNIEnv *env, jclass clazz, jdouble x, jdouble y)
{
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
        return log(exp((x - base) * LOG_2) + exp((y - base) * LOG_2)) * ONE_LOG_2 + base;
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
        
        return fastlog2(tot) + base;
}

JNIEXPORT jdouble JNICALL Java_net_derkholm_nmica_maths_NativeMath_addLog2___3D
  (JNIEnv *env, jclass clazz, jdoubleArray x)
{
	jdouble *elements;
	jsize length;
	jdouble tot;

	elements = (*env)->GetDoubleArrayElements(
        env,
        x,
        NULL
    );
    length = (*env)->GetArrayLength(env, x);
    
    tot = multiAddLog2(elements, length);
    
    (*env)->ReleaseDoubleArrayElements(env, x, elements, JNI_ABORT);
    
    return tot;
}
