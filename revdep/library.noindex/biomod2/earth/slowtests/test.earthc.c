// test.c: main() for testing earth c routines
// Comments containing "TODO" mark known issues

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <crtdbg.h>
#include "../../src/earth.h"

#define PRINT_TIME 0

#define sq(x)   ((x) * (x))

#if _MSC_VER            // microsoft
    // disable warning: 'vsprintf': This function or variable may be unsafe
    #pragma warning(disable: 4996)
#endif

//-----------------------------------------------------------------------------
void error(const char *args, ...)
{
    char s[1000];
    va_list p;
    va_start(p, args);
    vsprintf(s, args, p);
    va_end(p);
    printf("\nError: %s\n", s);
    // printf("Forcing a crash for the JIT debugger\n");
    // fflush(stdout);
    // _sleep(1000); // wait for fflush to finish
    // static volatile int* p999 = 0;
    // *p999 = 999;
    exit(-1);
}

//-----------------------------------------------------------------------------
static double getGcv(const int nTerms, // nbr basis terms including intercept
                const int nCases, double Rss, const double penalty)
{
    double cost;
    if (penalty < 0)    // special case: terms and knots are free
        cost = 0;
    else
        cost = nTerms + (penalty * (double)(nTerms-1) / 2);   // nKnots = (nTerms-1)/2

    return Rss / (nCases * sq(1 - cost/nCases));
}

//-----------------------------------------------------------------------------
static double RandUniform(void) // uniform rand number from -1 to +1
{
    return (double)((rand() % 20000) - 10000) / 10000;
}

//-----------------------------------------------------------------------------
static double RandGauss(void)   // standard normal random number
{
    double r = 0;
    for (int i = 0; i < 12; i++)    // by central limit theorem sum of uniforms is gaussian
        r += RandUniform();
    return r / 2;
}

//-----------------------------------------------------------------------------
static double funcNoise(const double x[], const int iResponse) { return RandGauss(); }

static double func0(const double x[], const int iResponse) { return x[0]; }

static double func1clean(const double x[], const int iResponse) { return x[0] + x[1]; }

static double func1(const double x[], const int iResponse) { return x[0] + x[1] + .1 * RandGauss(); }

static double func2(const double x[], const int iResponse) { return x[0] + x[1] + x[0]*x[1]; }

static double func3(const double x[], const int iResponse) { return cos(x[0]) + x[1]; }

static double func4(const double x[], const int iResponse) { return sin(2 * x[0]) + 2*x[1] + 0.5*x[0]*x[1]; }

static double func5(const double x[], const int iResponse) { return x[0] + x[1] + x[3] + x[4] + x[5] + x[1]*x[2] + (x[3]+1)*(x[4]+1)*x[5]; }

static double func4lin(const double x[], const int iResponse) { return x[0] + x[1] + x[3] + x[4]; }

static double func6(const double x[], const int iResponse) {                 // 5 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5]
            + .1 * RandGauss();
}

static double func6clean(const double x[], const int iResponse) {            // 5 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5];
}

static double func7(const double x[], const int iResponse) {                 // 10 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +x[6] +x[7] +x[8] +x[9] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5] + x[6]*x[7] + x[8]*x[9];
}

static double func8(const double x[], const int iResponse) {                 // 20 preds, 2nd order
    return  x[0] +x[1]+ x[2] +x[3] +x[4] +x[5] +x[6] +x[7] +x[8] +x[9] +
           x[10]+x[11]+x[12]+x[13]+x[14]+x[15]+x[16]+x[17]+x[18]+x[19] +
            x[0]*x[1] + x[2]*x[3] + x[4]*x[5] + x[6]*x[7] + x[8]*x[9] +
            + .1 * RandGauss();
}

static double func9(const double x[], const int iResponse) { return x[1]; }

static double func56(const double x[], const int iResponse) {    // Friedman MARS paper eqn 56
    return 0.1 * exp(4*x[0]) + 4 / (1 + exp(-20*(x[1]-0.5)) + 3*x[2] + 2*x[3] + x[4] + RandGauss());
}

// functions for testing multiple responses

static double func0_1clean(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func0(x, iResponse);
    else if (iResponse == 1)
        return func1clean(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func0_1(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func0(x, iResponse);
    else if (iResponse == 1)
        return func1(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func2_2(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func2(x, iResponse);
    else if (iResponse == 1)
        return func2(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func0_4(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func0(x, iResponse);
    else  if (iResponse == 1)
        return func4(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func0_2_4(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func0(x, iResponse);
    else  if (iResponse == 1)
        return func2(x, iResponse);
    else  if (iResponse == 2)
        return func4(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func2_4_0(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func2(x, iResponse);
    else  if (iResponse == 1)
        return func4(x, iResponse);
    else  if (iResponse == 2)
        return func0(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func4_2_0(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func4(x, iResponse);
    else  if (iResponse == 1)
        return func2(x, iResponse);
    else  if (iResponse == 2)
        return func0(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

static double func4_6(const double x[], const int iResponse)
{
    if (iResponse == 0)
        return func4(x, iResponse);
    else  if (iResponse == 1)
        return func6(x, iResponse);
    else
        error("bad iResponse");
    return 0;
}

// functions for testing NewVarPenalty

static double func1collinear(const double x[], const int iResponse)
    { return x[0] + x[1] + .001 * RandGauss(); }

static double func2collinear(const double x[], const int iResponse)
    { return cos(x[0]) + cos(x[1]) + .1 * RandGauss(); }

int cmp(const void *x, const void *y) // for qsort
{
    double x0 = *(double*)x;
    double y0 = *(double*)y;
    if (x0 < y0) return -1;
    if (x0 > y0) return  1;
    return 0;
}

//-----------------------------------------------------------------------------
static void TestEarth(char sTestName[],
                double (__cdecl *func)(const double xrow[], const int iResponse),
                const int nCases, const int nResponses, const int nPreds,
                const int nMaxDegree, const int nMaxTerms,
                const double Trace, const bool Format,
                const double ForwardStepThresh,
                const int nFastK, const double FastBeta, const double NewVarPenalty,
                const int seed,
                const double Collinear = 0) // used for testing NewVarPenalty
{
    #define y_(i,iResponse) y[(i) + (iResponse)*(nCases)]

    int *LinPreds  = (int *)calloc(nPreds, sizeof(int));

    double *x         = (double *)malloc(nCases    * nPreds     * sizeof(double));
    double *y         = (double *)malloc(nCases    * nResponses * sizeof(double));
    double *bx        = (double *)malloc(nCases    * nMaxTerms  * sizeof(double));
    bool   *BestSet   = (bool *)  malloc(nMaxTerms *              sizeof(bool));
    int    *Dirs      = (int *)   malloc(nMaxTerms * nPreds     * sizeof(int));
    double *Cuts      = (double *)malloc(nMaxTerms * nPreds     * sizeof(double));
    double *Residuals = (double *)malloc(nCases    * nResponses * sizeof(double));
    double *Betas     = (double *)malloc(nMaxTerms * nResponses * sizeof(double));

    static int nTest;
    nTest++;

    printf("=============================================================================\n");
    printf("TEST %d: %s n=%d p=%d\n", nTest, sTestName, nCases, nPreds);

    // init x
    srand(seed);
    int i;
    for (i = 0; i < nCases; i++)
        for (int iPred = 0; iPred < nPreds; iPred++) {
            double xtemp;
            xtemp = (double)((rand() % 20000) - 10000) / 10000;    // rand number from -1 to +1
            x[i + iPred * nCases] = xtemp;
        }

    // sort the first column of x, makes debugging easier
    qsort(x, nCases, sizeof(double), cmp);

    if (Collinear > 0) {
        // copy column 0 to 1 with added noise
        for (i = 0; i < nCases; i++)
            x[i + 1 * nCases] = x[i] + Collinear * RandGauss();
    }

    // init y
    double *xrow = (double *)malloc(nPreds * sizeof(double));
    for (i = 0; i < nCases; i++) {
        for (int iPred = 0; iPred < nPreds; iPred++)
            xrow[iPred] = x[i + iPred * nCases];
        for (int iResponse = 0; iResponse < nResponses; iResponse++)
            y_(i, iResponse) = func(xrow, iResponse);
    }
    free(xrow);

    double BestGcv;
    int nTerms, iReason;
    const double Penalty = ((nMaxDegree>1)? 3:2);
    clock_t Time = clock();
    if(Trace >= 4) {
        if(nResponses != 1)
            error("cannot use Trace>=4 with nResponses!=1");
        printf("           y");
        for(int iPred = 0; iPred < nPreds; iPred++)
            printf("         x%d", iPred);
        printf("\n");
        for(int i = 0; i < nCases; i++) {
            printf("%4d % 7.5f", i, y[i]);
            for(int iPred = 0; iPred < nPreds; iPred++) {
                printf(" % 7.5f", x[i + iPred * nCases]);
            }
            printf("\n");
        }
        printf("\n");
    }
    Earth(&BestGcv, &nTerms, &iReason, BestSet, bx, Dirs, Cuts, Residuals, Betas,
        x, y, NULL, // weights are NULL
        nCases, nResponses, nPreds, nMaxDegree, nMaxTerms, Penalty, ForwardStepThresh,
        0, 0,   // MinSpan, EndSpan
        true,   // Prune
        nFastK, FastBeta, NewVarPenalty, LinPreds,
        2 /*AdjustEndSpan*/, true /*AutoLinPred*/, true /*UseBetaCache*/,
        Trace, NULL);

    // calc nUsedTerms

    int nUsedTerms = 0;
    for (int iTerm = 0; iTerm < nTerms; iTerm++)
        if (BestSet[iTerm])
            nUsedTerms++;

    // calc RSquared, GRSquared

    for (int iResponse = 0; iResponse < nResponses; iResponse++) {
        double Rss = 0, Tss = 0, meanY = 0;
        for (i = 0; i < nCases; i++)
            meanY += y_(i, iResponse);
        meanY /= nCases;
        xrow = (double *)malloc(nPreds * sizeof(double));
        double *yHat = (double *)malloc(nResponses * sizeof(double));
        for (i = 0; i < nCases; i++) {
            for (int iPred = 0; iPred < nPreds; iPred++)
                xrow[iPred] = x[i + iPred * nCases];
            PredictEarth(yHat, xrow, BestSet, Dirs, Cuts, Betas, nPreds, nResponses, nTerms, nMaxTerms);
            double Residual = y_(i,iResponse) - yHat[iResponse];
            Rss += sq(Residual);
            Tss += sq(y_(i,iResponse) - meanY);
        }
        free(yHat);
        free(xrow);
        const double RSq =  1 - Rss/Tss;
        const double GcvNull =  getGcv(1, nCases, Tss, Penalty);
        const double GRSq =  1 - getGcv(nUsedTerms, nCases, Rss, Penalty) / GcvNull;

#if PRINT_TIME
        double TimeDelta = (double)(clock() - Time) / CLOCKS_PER_SEC;
#else
        double TimeDelta = 99.99;
#endif
        // show results
        if (nResponses > 1) {
            printf("RESULT %d Response %d: GRSq %.5g RSq %.5g nTerms %d of %d of %d",
                   nTest, iResponse+1, GRSq, RSq,
                   nUsedTerms, nTerms, nMaxTerms);
            if (iResponse == 0)
                printf(" FUNCTION %s n=%d p=%d [%.2f secs]",
                   sTestName, nCases, nPreds, TimeDelta);
            printf("\n");
        }
        else
            printf("RESULT %d: GRSq %g RSq %g nTerms %d of %d of %d "
                   "FUNCTION %s n=%d p=%d [%.2f secs]\n",
                   nTest, GRSq, RSq, nUsedTerms, nTerms, nMaxTerms,
                   sTestName, nCases, nPreds, TimeDelta);
    }
    if (Format && Trace != 0) {
        printf("\nTEST %d: FUNCTION %s n=%d p=%d\n", nTest, sTestName, nCases, nPreds);
        FormatEarth(BestSet, Dirs, Cuts, Betas, nPreds, nResponses, nTerms, nMaxTerms, 3, 1e-6);
        printf("\n");
    }
    free(LinPreds);
    free(x);
    free(y);
    free(BestSet);
    free(Dirs);
    free(Cuts);
    free(Residuals);
    free(Betas);
    free(bx);
}

//-----------------------------------------------------------------------------
int main(void)
{
  clock_t Time = clock();
                                       // func  nCases     nResp nPreds  nMaxDegree nMaxTerms  Trace Form Thresh K B N s
  TestEarth("noise",                   funcNoise, 1000,        1,     1,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("x0",                      func0,       10,        1,     1,          2,       51,     7,true,0.001,20,1,0,99);

  // intercept only models
  TestEarth("x0",                      func0,       10,        1,     1,          2,        1,     3,true,0.001,20,1,0,99);
  TestEarth("x0",                      func0,       10,        1,     1,          2,        2,     3,true,0.001,20,1,0,99);

  TestEarth("x0",                      func0,     1000,        1,     1,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("x0 + noise",              func0,     1000,        1,   1+1,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("x0 + x1",                 func1,     1000,        1,     2,          2,       11,     7,true,0.001,20,1,0,99);
  TestEarth("x0 + x1 + noise",         func1,     1000,        1,   2+8,          2,       51,     0,true,0.001,20,1,0,99);
  TestEarth("x0 + x1 + x0*x1",         func2,       30,        1,     2,          2,       51,     4,true,0.001,20,1,0,99);
  TestEarth("x0 + x1 + x0*x1",         func2,     1000,        1,     2,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("x0 + x1 + x0*x1",         func2,     1000,        1,     2,          2,       51,   1.5,true,0.001,20,1,0,99);
  TestEarth("cos(x0) + x1",            func3,     1000,        1,     2,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("sin(2*x0)+2*x1*.5*x0*x1", func4,     1000,        1,     2,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("sin(2*x0)+2*x1*.5*x0*x1", func4,     1000,        1,     3,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("3rd order, mi=2 ni=11",   func5,     1000,        1,     6,          2,       11,     1,true,0.001,20,1,0,99);
  TestEarth("3rd order, mi=2 ni=51",   func5,     1000,        1,     6,          2,       51,     2,true,0.001,20,1,0,99);
  TestEarth("3rd order, mi=3",         func5,     1000,        1,     6,          3,       51,     3,true,0.001,20,1,0,99);
  TestEarth("5 preds + noise",         func6,      200,        1,  5+10,          2,      101,     3,true,0.001,20,1,0,99);
  TestEarth("5 preds clean",           func6clean, 200,        1,  5+10,          2,      101,     3,true,0.001,20,1,0,99);
  TestEarth("10 preds + noise",        func7,      200,        1, 10+40,          2,      101,     3,true,0.001,20,1,0,99);
  TestEarth("20 preds + noise,",       func8,      100,        1, 20+10,          2,      101,     3,true,0.001,20,1,0,99);
  TestEarth("20 preds + noise,",       func8,      400,        1, 20+10,          2,      101,     3,true,0.001,20,1,0,99);
  TestEarth("3rd order, mi=3 + noise", func5,     1000,        1,    10,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("eqn56 mi=1",              func56,     300,        1,     6,          1,      101,     3,true,0.001,20,1,0,99);
  TestEarth("eqn56 mi=2",              func56,     300,        1,     6,          2,       51,     3,true,0.001,20,1,0,99);
  TestEarth("eqn56 mi=10",             func56,     300,        1,     6,         10,       51,     3,true,0.001,20,1,0,99);
  // Following two tests are slow so are commented out (take more than half a second each)
  // TestEarth("eqn56 mi=10",          func56,    1000,        1,     6,         10,      101,     3,true,0.001,20,1,0,99);
  // TestEarth("eqn56 mi=10",          func56,    5000,        1,     6,         10,      101,     3,true,0.001,20,1,0,99);

  TestEarth("x0 + x1 + x0*x1",         func2,       30,        1,     2,          2,       51,     3,true,0.001,99,1,0,99);
  TestEarth("x0 + x1 + x0*x1",         func2,       30,        1,     2,          2,       51,     3,true,0.001, 4,0,0,99);
  TestEarth("x0 + x1 + x0*x1",         func2,       30,        1,     2,          2,       51,     3,true,0.001, 4,1,0,99);

  // test multiple responses                func    nCases     nResp nPreds  nMaxDegree nMaxTerms  Trace Form  Thresh K B N s

  TestEarth("x0|x0+x1 degree=1",     func0_1clean,  30,        2,     2,          1,       51,     3, true,0.001,20,1,0,99);

  TestEarth("x0|x+x1+noise",         func0_1,      100,        2,     2,          1,       51,     3, true,0.001,20,1,0,99);

  TestEarth("x0+x1+x0*x1|x0+x1+x0*x1 degree=1",
                                     func2_2,      100,        2,     2,          1,       51,     3, true,0.001,20,1,0,99);

  TestEarth("x0+x1+x0*x1|x0+x1+x0*x1 degree=2",
                                     func2_2,      100,        2,     2,          2,       51,     3, true,0.001,20,1,0,99);

  TestEarth("x0|sin(2*x0) + 2*x1 + 0.5*x0*x1 + 8 noise preds",
                                     func0_4,      200,        2,    10,          2,      101,     3, true,0.001,20,1,0,99);

  TestEarth("x0|x0+x1+x0*x1|sin(2*x0) + 2*x1 + 0.5*x0*x1  + 8 noise preds",
                                     func0_2_4,    200,        3,   3+8,          2,      101,     3, true,0.001,20,1,0,99);

  TestEarth("|x0+x1+x0*x1|sin(2*x0) + 2*x1 + 0.5*x0*x1|x0  + 8 noise preds",
                                     func2_4_0,    200,        3,   3+8,          2,      101,     3, true,0.001,20,1,0,99);

  TestEarth("sin(2*x0) + 2*x1 + 0.5*x0*x1|x0+x1+x0*x1|x0  + 8 noise preds",
                                     func4_2_0,    200,        3,   3+8,          2,      101,     3, true,0.001,20,1,0,99);

  //TODO following gives lousy GRSq for Response 2, investigate
  TestEarth("sin(2*x0) + 2*x1 + 0.5*x0*x1|2nd order 6 preds + noise",
                                     func4_6,     1000,        2,     6,         2,      101,      3, true,0.001,20,1,0,99);

  // test NewVarPenalty                     func    nCases     nResp nPreds  nMaxDegree nMaxTerms   Trace Form Thresh  K B  NP  s Colin

  // basis functions (after pruning) include both predictors
  TestEarth("cos(x1) + cos(x2), x1 and x2 xcollinear, NewVarPenalty=0",
                               func2collinear,      100,       1,     2,          1,       51,     3, true, 0.001,20,1,  0, 99,.07);

  // basis functions (after pruning) include only one predictor
  TestEarth("cos(x1) + cos(x2), x1 and x2 xcollinear, NewVarPenalty=0",
                               func2collinear,      100,       1,     2,          1,       51,     3, true, 0.001,20,1,0.1, 99,.07);

#if PRINT_TIME
  printf("[Total time %.2f secs]\n", (double)(clock() - Time) / CLOCKS_PER_SEC);
#endif

  return 0;
}
