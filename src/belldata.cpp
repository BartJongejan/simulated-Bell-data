/* belldata.cpp
Bart Jongejan, September 2025, Copenhagen, Denmark

Command line program to create
1) Data that could have been created in an ideal Bohm-Bell experiment
2) Hidden data that could have caused 1)

The hidden data in 2) is shared by Alice and Bob

command line:

    ./belldata <dimensions> <# Alice's keys> <# Bob's keys> <# trials> [<offset>]
e.g.
    ./belldata 3 30 30 30000000

An offset > 0 causes a mismatch between Alice's and Bob's observations. All statistical
correlations are absent if offset > 0.

The mentioned "keys" are the keys on an apparatus. Alice and Bob have each an apparatus.
Alice and Bob press keys with labels a1,a2, etc. and  b1,b2, etc.
A yellow lamp lights up when Alice and Bob have to check the colour {not-green,green}
that is visible in a slit on the apparatus.

The keys take the place of Bell's vectors, and the colours take the place of Bell's spin-up/down.

There are many more keys than there are vectors in Bell's original experiment.

In contrast to Bell's vectors, which are embedded in a reference system that Alice and Bob share,
the apparatuses are prepared without sharing anything.

With offset=0, the generated data houses a strong correlation between keys and colour observations.
Expectation values P(a,b) can be found that, once filled in into the CHSH inequality, can violate it.
For dimension=2, no violation occurs.
For higher dimensions, violations occur. The higher the dimension, the larger the violation, with a
ceiling of 4.

The correlations can be used (not demonstrated here) to construct the configuration of spatial
directions in Alice's and Bob's apparatuses. Thus, in the end, the spatial system of reference
that we ignore to begin with is regained.

This program works backwards: from the spatial configuration towards the observed data.
(1) It starts with determining two sets of uniform distributed vectors. One for Alice's apparatus,
    the other for Bob's.
(2) It determines the separation H(a,b) between each pair of vectors (a,b).
    The separation is a function 0 <= H <= 2 that has a uniform distribution if a and b are uniform
    distributed. The exact form depends on the number of dimensions:
    In 2 dimensions, H(a,b) is equal to the angle between a and b divided by π/2
    In 3 dimensions, H(a,b) is equal to 1-cos(angle(a,b))
    etc.
(3) For each trial, the program selects two keys a and b and also three random numbers a, b,
    and λ, each between 0 and 4, under the restriction that the separation between a and b is equal
    to H(a, b).
(4) The values H(a,λ) and H(b,λ) are computed. These are the hidden variables that decide
    the colour observation. A value >=1 is the colour green, a value < 1 is not-green.

The 'reference system' one can construct from the data is not Newton's absolute space.
The separation H takes two keys that are not local to each other. It cannot take two keys from the
same location. It is a complete bipartite graph rather than a full-fledged system of reference.
*/

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>

#define MONTECARLO 1 // 0: λ increases with fixed increments; 1: random 0 <= λ < 4 
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923

typedef int         trial;
typedef int         dimension;
typedef int         key_choice;
typedef int         binaryoutcome;
typedef int         count;
typedef double      betweenminus1and1;
typedef double      between0and2;
typedef double      angle;
typedef double      between0and4;
typedef double      CHSHsum;
typedef double      component;

typedef struct
    {
    key_choice      a;          // Alice's key choice
    key_choice      b;          // Bob's key choice
    trial           n;          // Number of experiments involving these key choices
    trial           naup;       // Number of Alice's positive results. Should be close to n/2
    trial           nbup;       // Number of Bob's positive results. Should be close to n/2
    trial           nabup;      // Number of experiments where both Alice and Bob obtain positive results.
    trial           nabdown;    // Number of experiments where both Alice and Bob obtain negative results. Should be close to nabup
    between0and2    Hab;        // In 3 Dimensions: 1 minus the inner product of Alice and Bob's SG instrument unit vectors
    between0and2    difference; // Should be close to 1-Hab
    } matrixcell;

typedef struct
    {
    trial           e;          // Initially a random integer. Later repurposed as trial number.
    key_choice      a;          // Alice's key choice
    key_choice      b;          // Bob's key choice
    between0and4    al;         // 0 <= aλ < 4, aλ = |λ-a| is the hidden variable, relative to Alice.
    between0and4    bl;         // 0 <= bλ < 4, bλ = |λ-b| is the hidden variable, relative to Bob.
    } result;

typedef struct
    {
    key_choice      NAliceKeys;
    key_choice      NBobKeys;
    matrixcell* cell;
    count           uu;
    count           dd;
    count           du;
    count           ud;
    } matrix;

static dimension    Dimensions = 4; // Set by first parameter on the command line
static FILE* fpstat;

class vectr
    {
    public:
    public:
        component* e;
        void initRandomUnitVector()
            {// https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
            double sum2 = 0.0;
            for(dimension k = 0; k < Dimensions; ++k)
                {
                int iU = 1 + rand();
                int iV = 1 + rand();
                double U = (double)(iU) / ((double)RAND_MAX + 1.0);
                double V = (double)(iV) / ((double)RAND_MAX + 1.0);
                double root = sqrt(-2.0 * log(U));
                component XXd = root * cos(2.0 * M_PI * V);
                e[k] = XXd;
                sum2 += XXd * XXd;
                }
            double f = sqrt(sum2);
            if(f > 0.0)
                {
                for(dimension col = 0; col < Dimensions; ++col)
                    e[col] /= f;
                }
            }
        void init()
            {
            dimension dim = Dimensions;
            if(dim < 3)
                dim = 3;
            e = new component[dim];
            for(dimension i = 0; i < dim; ++i)
                e[i] = 0.0;
            }
        vectr()
            {
            if(Dimensions > 0)
                init();
            else
                e = 0;
            };
        ~vectr()
            {
            delete[] e;
            }
        betweenminus1and1 dotProduct(vectr* b)
            {
            assert(b != NULL);
            betweenminus1and1 ret = 0;
            for(dimension c = 0; c < Dimensions; ++c)
                {
                ret += e[c] * b->e[c];
                }
            return ret;
            }
    };

class apparatus
    {
    public:
        char* name;     // e.g. Alice, Bob
        vectr* vect;
        key_choice N;   // Number of keys. e.g., 30
        apparatus(key_choice N, const char* name) :N(N)
            {
            key_choice i = 0;
            this->name = new char[strlen(name) + 1];
            strcpy(this->name, name);
            vect = new vectr[N];
            for(i = 0; i < N; ++i)
                {
                vect[i].initRandomUnitVector();
                }
            }
        ~apparatus()
            {
            delete[] name;
            delete[] vect;
            }
    };

struct StandardDev
    {
    double mean;
    double E;
    double E2;
    double min, max;
    count N;
    StandardDev() :mean(-1), E(0), E2(0), min(DBL_MAX), max(-DBL_MAX), N(0)
        {}
    void add(double val, int pass)  // Must be called twice. Mean is computed after first pass.
                                    // In second pass the square of E is computed.
        {
        if(pass == 2)
            {
            if(mean == -1)
                mean = E / (double)N;
            val -= mean;
            E2 += val * val;
            }
        else
            {
            E += val;
            if(val < min) min = val;
            if(val > max) max = val;
            ++N;
            }
        }
    count n() { return N; }
    double avg()
        {
        return E / (double)N;
        }
    double stddev(FILE* fpstat)
        {
        double SD = -1;
        if(N > 0)
            {
            double S = E2 / (double)N;
            printf("min %f, max %f, N %d, E %f, E2 %f,  E2/N %f \n", min, max, N, E, E2, S);
            if(S > 0)
                SD = sqrt(S);
            else if(S == 0)
                SD = 0;
            else
                SD = -1.0;
            printf("SD %f\n\n", SD);
            if(fpstat)
                fprintf(fpstat, "%f", SD);
            }
        return SD;
        }
    void print(FILE* fp, const char* txt)
        {
        fprintf(fp, "%f +/- ", avg());
        stddev(fp);
        fprintf(fp, "   N %d   min %f  median %f  max %f   %s\n", N, min, (max - min) / 2.0, max, txt);
        }
    };

static void initMatrix(matrix* Matrix, key_choice SizeOfAlicesSampleSpace, key_choice SizeOfBobsSampleSpace)
    {
    Matrix->dd = Matrix->uu = Matrix->du = Matrix->ud = 0; // These are counters, so they start at 0.
    Matrix->NAliceKeys = SizeOfAlicesSampleSpace;
    Matrix->NBobKeys = SizeOfBobsSampleSpace;
    Matrix->cell = new matrixcell[SizeOfAlicesSampleSpace * SizeOfBobsSampleSpace];
    }

static betweenminus1and1 dotProduct(apparatus* AliceApparatus, apparatus* BobApparatus, key_choice p, key_choice q, dimension DIM)
    {
    return AliceApparatus->vect[p].dotProduct(BobApparatus->vect + q);
    }

/* :
Compute hyperspherical cap using a recursive algorithm.
*/
static angle I(angle gamma, dimension n)
    {
    angle result = 0;
    if(n == 0)
        result = gamma;
    else if(n == 1)
        result = 1.0 - cos(gamma);
    else
        {
        double dn = (double)n;
        result = -1.0 / dn * pow(sin(gamma), dn - 1.0) * cos(gamma) + (dn - 1.0) / dn * I(gamma, n - 2);
        }
    return result;
    }

static between0and2 E(angle gamma, dimension n)
    {
    static angle PI = M_PI;
    between0and2 result;
    result = 2.0 * I(gamma, n) / I(PI, n);
    return result;
    }

static between0and2 separation(betweenminus1and1 cosin, dimension n)
    {
    angle gamma = acos(cosin);
    return E(gamma, n);
    }

static matrixcell* theCell(matrix Matrix, key_choice a, key_choice b)
    {
    return Matrix.cell + b * Matrix.NAliceKeys + a;
    }

static void printpair(FILE* fp, matrix Matrix, result* Pr, apparatus* NAliceKeys, apparatus* NBobKeys)
    {
    between0and2 aFluid = Pr->al > 2.0 ? 4.0 - Pr->al : Pr->al;
    between0and2 bFluid = Pr->bl > 2.0 ? 4.0 - Pr->bl : Pr->bl;
    fprintf(fp, "%d\ta%d\tb%d\t%f\t%f\t%f\t%f\t%f\t%s\t%s\n",
            Pr->e,
            1 + Pr->a, // 1+ because lowest key label is a1, not a0
            1 + Pr->b, // 1+ because lowest key label is b1, not b0
            theCell(Matrix, Pr->a, Pr->b)->Hab,
            Pr->al, Pr->bl,
            aFluid, bFluid,
            aFluid < 1.0 ? "black" : "green", bFluid < 1.0 ? "black" : "green");
    }

static void printpairNarrow(FILE* fp, matrix Matrix, result* Pr, apparatus* NAliceKeys, apparatus* NBobKeys)
    { /* Compact output, but many lines and therefore very large all the same. */
    between0and2 aFluid = Pr->al > 2.0 ? 4.0 - Pr->al : Pr->al;
    between0and2 bFluid = Pr->bl > 2.0 ? 4.0 - Pr->bl : Pr->bl;
    key_choice A = Pr->a;
    key_choice B = Pr->b;
    char a = A < 10 ? '0' + A : 'A' - 10 + A;
    char b = B < 10 ? '0' + B : 'A' - 10 + B;
    fprintf(fp, "%c%c%c%c\n", a, b, aFluid < 1.0 ? '-' : '+', bFluid < 1.0 ? '-' : '+');
    }

static int randomizeLambdas(const void* x, const void* y)
    {
    trial X = (*(result**)x)->e;
    trial Y = (*(result**)y)->e;
    if(X > Y)
        return 1;
    else if(X < Y)
        return -1;
    else
        return 0;
    }

static bool gt(CHSHsum a, CHSHsum b)
    {
    return a > b;
    }

static bool lt(CHSHsum a, CHSHsum b)
    {
    return a < b;
    }

/* Compute normalized correlations that are in the interval [0,2]
   0.0 : perfect correlation
   1.0 : no correlation
   2.0 : perfect anti-correlation
*/
static void normalizeCorrelations(matrix Matrix)
    {
    key_choice p;
    key_choice q;
    for(p = 0; p < Matrix.NAliceKeys; ++p)
        for(q = 0; q < Matrix.NBobKeys; ++q)
            {
            matrixcell* Cell = theCell(Matrix, p, q);
            if(Cell->n > 0)
                {
                Cell->difference = 2.0 * (double)(Cell->n - Cell->nabup - Cell->nabdown) / (double)Cell->n;
                }
            else
                Cell->difference = 0;
            }
    }

static void createLambdas(result* Experiment, trial Ntrials)
    {
    /*
    lambda climbs monotonously from 0 to 4 in Ntrials small, equal steps in the pLambda array.
    The e field is used to disorganise the experiments, so that lambda varies randomly from
    trial to the next trial in the Experiment array. This is not necessary, but adds to the
    photo-realism of the produced 'data'.
    */
    trial t;
    result** pLambda = new result * [Ntrials];
    for(t = 0; t < Ntrials; ++t)
        {
        pLambda[t] = Experiment + t;
        pLambda[t]->e = rand();
        }
    /* Disorganise lambda by sorting the column e that contains quasi-random numbers. */
    qsort(pLambda, Ntrials, sizeof(result*), randomizeLambdas);

    /* Set lowest value for λ. This is not the initial value of λ! The array of λ is randomized.*/
    between0and4 lambda = 0;
#if !MONTECARLO
    between0and4 deltalambda = 4.0 / (double)Ntrials;
    lambda = 0.5 * deltalambda;
#endif

    /* The uniform distribution of λ can be continuous or discrete. */
    for(trial L = 0; L < Ntrials; ++L) // Run over Experiments in random fashion.
        {
        pLambda[L]->al = lambda - 0; /* a = 0 */
#if MONTECARLO
        /* Continuous uniform distribution of λ */
        lambda = (between0and4)rand() * (4.0 / (double)RAND_MAX);
#else
        /* Discrete uniform distribution of λ */
        /* Note that the λ come in randomized order over the trials! */
        lambda += deltalambda;
#endif
        }
    delete[] pLambda;
    }

static void createChoices(key_choice* Choice, key_choice SizeOfSampleSpace, trial Ntrials)
    {
    /* Determine the keys Alice and Bob press in each trial of the experiment */
    trial t;
    for(t = 0; t < Ntrials; ++t)
        Choice[t] = rand() % SizeOfSampleSpace;
    }

static void determineHab(matrix CorrelationMatrix, dimension dim, apparatus* AliceApparatus, apparatus* BobApparatus)
    {
    key_choice p;
    key_choice q;
    matrixcell* Cell;
    key_choice SizeOfAlicesSampleSpace = AliceApparatus->N;
    key_choice SizeOfBobsSampleSpace = BobApparatus->N;
    for(p = 0; p < SizeOfAlicesSampleSpace; ++p)
        {
        for(q = 0; q < SizeOfBobsSampleSpace; ++q)
            {
            betweenminus1and1 in;
            in = dotProduct(AliceApparatus, BobApparatus, p, q, dim);
            Cell = theCell(CorrelationMatrix, p, q);
            Cell->a = p;
            Cell->b = q;
            Cell->n = 0;
            Cell->nabdown = 0;
            Cell->nabup = 0;
            Cell->naup = 0;
            Cell->nbup = 0;
            Cell->Hab = separation(in, dim - 2); /* The most important line! */
            Cell->difference = 0.0;
            }
        }
    }

static void doTrials(result* Experiment, trial Ntrials, key_choice* AlicesChoice, key_choice* BobsChoice, matrix CorrelationMatrix)
    {
/*
* Create the table for trials (pairs of measurements).
* Each of the Ntrials records contains two key choices and two lambdas,
* one for each key choice. The records are numbered 0 ... (Ntrials - 1)
*/
    trial t;
    for(t = 0; t < Ntrials; ++t)
        {
        between0and4 al;
        between0and4 bl;
        between0and2 b;
        between0and2 Hab; // separation between a and b
        result* theTrial = Experiment + t;
        key_choice AlicesKey = AlicesChoice[t];
        key_choice BobsKey = BobsChoice[t];
        matrixcell* Cell = theCell(CorrelationMatrix, AlicesKey, BobsKey);
        al = theTrial->al;  // 0 < aλ < 4, aλ = |λ-a|, a=0, 0<= λ < 4 is a stochastic variable
        Hab = Cell->Hab;    // In 3 Dimensions, Hab = H(a,b) is 1-(inner product of two unit vectors).
        /* We set a=0 and b = H(a,b) = separation between a and b,
           where H(a,b) is uniformly distributed in the interval [0,2] over all trials. */
        b = Hab;            // 0 <= b <= 2

        bl = al - b;
        if(bl < 0.0)
            bl += 4.0;
        else if(bl > 4.0)
            bl -= 4.0;
        theTrial->bl = bl; // 0 <= bλ < 4, bλ = |λ - b|, 0 <= b <= 2, 0<= λ < 4 is a stochastic variable
        theTrial->e = t;   // Overwrites random number that already has been used to randomize λ.
        theTrial->a = AlicesKey;
        theTrial->b = BobsKey;
        }
    }

static void count_SpinUp_SpinDown(matrix* CorrelationMatrix, result* Experiment, trial Ntrials, trial noff)
    {
    for(trial t = 0; t < Ntrials; ++t)
        {
        between0and4 al;  // λ
        between0and4 bl;
        between0and2 diffah;  // |a-λ|
        between0and2 diffbh;  // |b-λ|
        result* AlicesExperiment = Experiment + t;
        result* BobsExperiment = Experiment + ((t + noff) % Ntrials);
        key_choice AlicesSetting = AlicesExperiment->a;
        key_choice BobsSetting = BobsExperiment->b;
        al = AlicesExperiment->al;  // 0 < λ < 4, λ is a stochastic variable
        bl = BobsExperiment->bl;  // 0 < λ < 4, λ is a stochastic variable

        if(al > 2.0)
            diffah = 4.0 - al;
        else
            diffah = al;

        if(bl > 2.0)
            diffbh = 4.0 - bl;
        else
            diffbh = bl;

        matrixcell* Cell = theCell(*CorrelationMatrix, AlicesSetting, BobsSetting);
        if(diffah < 1.0 && diffbh < 1.0) ++CorrelationMatrix->uu;
        if(diffah < 1.0 && diffbh >= 1.0) ++CorrelationMatrix->ud;
        if(diffah >= 1.0 && diffbh < 1.0) ++CorrelationMatrix->du;
        if(diffah >= 1.0 && diffbh >= 1.0) ++CorrelationMatrix->dd;
        if(diffah < 1.0)
            {
            ++(Cell->naup);
            if(diffbh < 1.0)
                ++(Cell->nabup); // both up
            }
        else
            {
            if(diffbh > 1.0)
                {
                ++(Cell->nabdown); // both down
                }
            }

        if(diffbh < 1.0)
            ++(Cell->nbup);

        ++(Cell->n);
        }
    printf("Alice green %d Bob green %d all %d\n"
           , CorrelationMatrix->du + CorrelationMatrix->dd
           , CorrelationMatrix->ud + CorrelationMatrix->dd
           , CorrelationMatrix->du + CorrelationMatrix->uu + CorrelationMatrix->ud + CorrelationMatrix->dd);
    }

static CHSHsum CHSHtheoreticalExtreme(dimension DIM)
    {
    static angle alpha = M_PI / 4.0;
    angle tripleAlpha;
    tripleAlpha = 3.0 * alpha;
    CHSHsum chsh;
    chsh = 3.0 * E(alpha, DIM - 2) - E(tripleAlpha, DIM - 2);
    fprintf(fpstat, "Theoretical maximum under/overshooting, based on optimal angle (not actual settings of Alice and Bob) %f\n", chsh);
    return chsh;
    }

static CHSHsum CHSHexperimetalExtreme(key_choice SizeOfAlicesSampleSpace, key_choice SizeOfBobsSampleSpace, matrix Matrix,
                                      bool (*comp)(CHSHsum a, CHSHsum b))
    {
    key_choice A0 = 0, A1 = 0, B0 = 0, B1 = 0;
    CHSHsum max;
    between0and2 Cors[4];
    key_choice a0 = 0, a1 = 0, b0 = 0, b1 = 0;
    int g;
    if(comp == gt)
        max = -10.0;
    else
        max = 10.0;
    int min = -1;
    for(a1 = 1; a1 < SizeOfAlicesSampleSpace; ++a1)
        {
        for(b1 = 1; b1 < SizeOfBobsSampleSpace; ++b1)
            {
            for(a0 = 0; a0 < a1; ++a0)
                {
                for(b0 = 0; b0 < b1; ++b0)
                    if((theCell(Matrix, a0, b0)->n > 3)
                       && (theCell(Matrix, a1, b0)->n > 3)
                       && (theCell(Matrix, a1, b1)->n > 3)
                       && (theCell(Matrix, a0, b1)->n > 3))
                        {
                        between0and2 S;
                        between0and2 cors[4];
                        cors[0] = theCell(Matrix, a0, b0)->difference;
                        cors[1] = theCell(Matrix, a1, b0)->difference;
                        cors[2] = theCell(Matrix, a1, b1)->difference;
                        cors[3] = theCell(Matrix, a0, b1)->difference;

                        S = cors[0];
                        S += cors[1];
                        S += cors[2];
                        S += cors[3];
                        for(g = 0; g < 4; ++g)
                            {
                            CHSHsum chsh;
                            chsh = (S - 2.0 * cors[g]);
                            if(comp(chsh, max))
                                {
                                max = chsh;
                                A1 = a1; A0 = a0; B1 = b1; B0 = b0;
                                for(int f = 0; f < 4; ++f)
                                    {
                                    Cors[f] = cors[f] - 1.0;
                                    if(f == g)
                                        {
                                        min = f;
                                        }
                                    }
                                }
                            }
                        }
                }
            }
        }
    double sum = 0.0;
    printf("a%d b%d: %f\n", A0 + 1, B0 + 1, Cors[0] + 1.0);
    printf("a%d b%d: %f\n", A1 + 1, B0 + 1, Cors[1] + 1.0);
    printf("a%d b%d: %f\n", A0 + 1, B1 + 1, Cors[2] + 1.0);
    printf("a%d b%d: %f\n", A1 + 1, B1 + 1, Cors[3] + 1.0);
    for(int f = 0; f < 4; ++f)
        {
        if(f == min)
            {
            sum -= Cors[f];
            printf("(%f) ", Cors[f] + 1);
            printf("-");
            }
        else
            {
            sum += Cors[f];
            printf("(%f) ", Cors[f] + 1);
            printf("+");
            }
        printf("%f ", Cors[f]);
        }
    printf("= (%f) %f \n", sum + 2.0, sum);
    fprintf(fpstat, "CHSH extreme, experimentally:%f a%d a%d b%d b%d\n",
            max, A0 + 1, A1 + 1, B0 + 1, B1 + 1);
    return max;
    }

static void writeResults(matrix CorrelationMatrix, apparatus* AliceApparatus, apparatus* BobApparatus, trial Ntrials, dimension dim, result* Experiment)
    {
    char name[256];
    sprintf(name, "HVsimulation-dim%d-Alice%d-Bob%d-iterations%d.tab", dim, AliceApparatus->N, BobApparatus->N, Ntrials);
    FILE* fp = fopen(name, "w");
    if(fp)
        {
        printf("short\n");
        fprintf(fp, "#\tAlice's key\tBob's key\tHab\taλ\tbλ\taFluid\tbFluid\tA\tB\n");
        trial rows = Ntrials;
        if(rows > 10000)
            rows = 10000;
        for(trial t = 0; t < rows; ++t)
            {
            printpair(fp, CorrelationMatrix, Experiment + t, AliceApparatus, BobApparatus);
            }
        fclose(fp);
        }
#if 0 /* Set to non-zero value if you need output from all trials. Output will be BIG! */
    sprintf(name, "HVsimulation-dim%d-Alice%d-Bob%d-iterations%d-BIG.tab", dim, AliceApparatus->N, BobApparatus->N, Ntrials);
    fp = fopen(name, "w");
    if(fp)
        {
        printf("long\n");
        for(trial t = 0; t < Ntrials; ++t)
            {
            printpairNarrow(fp, CorrelationMatrix, Experiment + t, AliceApparatus, BobApparatus);
            }
        fclose(fp);
        }
#endif
    }

static void doExperiment(dimension dim, key_choice NAlice, key_choice NBob, trial Ntrials, trial noff)
    {
    /* Create two apparatuses, each with its own set of random settings. */
    srand(0);
    apparatus* AliceApparatus = new apparatus(NAlice, "Alice");
    apparatus* BobApparatus = new apparatus(NBob, "Bob");
    srand(0);
    char name[256];
    sprintf(name, "statistics%dx%d-%d-%dD.txt", AliceApparatus->N, BobApparatus->N, Ntrials, dim);
    fpstat = fopen(name, "wb");
    key_choice* AlicesChoice = new key_choice[Ntrials];// AlicesChoice[e] is a stochastic variable. It represents the decision Alice makes in the trial with number e.
    key_choice* BobsChoice = new key_choice[Ntrials];//   BobsChoice[e] is a stochastic variable. It represents the decision Bob   makes in the trial with number e.
    result* Experiment = new result[Ntrials];
    /* Initilise a matrix with a cell for every combination of key choices by Alice and Bob */
    matrix CorrelationMatrix;
    initMatrix(&CorrelationMatrix, AliceApparatus->N, BobApparatus->N);
    createChoices(AlicesChoice, AliceApparatus->N, Ntrials);
    createChoices(BobsChoice, BobApparatus->N, Ntrials);
    createLambdas(Experiment, Ntrials);
    determineHab(CorrelationMatrix, dim, AliceApparatus, BobApparatus);
    doTrials(Experiment, Ntrials, AlicesChoice, BobsChoice, CorrelationMatrix);
    count_SpinUp_SpinDown(&CorrelationMatrix, Experiment, Ntrials, noff);
    normalizeCorrelations(CorrelationMatrix);

    fprintf(fpstat, "\nThe extrema of the CHSH expression for %d Dimensions:\n", dim);
    CHSHtheoreticalExtreme(dim);
    CHSHexperimetalExtreme(AliceApparatus->N, BobApparatus->N, CorrelationMatrix, gt);
    CHSHexperimetalExtreme(AliceApparatus->N, BobApparatus->N, CorrelationMatrix, lt);
    writeResults(CorrelationMatrix, AliceApparatus, BobApparatus, Ntrials, dim, Experiment);
    fclose(fpstat);
    delete[] AlicesChoice;
    delete[] BobsChoice;
    delete[] Experiment;
    delete[] CorrelationMatrix.cell;
    delete AliceApparatus;
    delete BobApparatus;
    }

int main(int argc, char* argv[])
    {
    printf("argc %d\n", argc);
    if(argc < 3)
        {
        printf("belldata <N-Dimensions> <N-Alice-settings> <N-Bob-settings> <N-experiments> [<N-Bob-offset>]\n(x resp. y, NAlice and NBob > 3)\n");
        return 1;
        }
    Dimensions = strtol(argv[1], NULL, 10);
    if(argc < 4)
        {
        printf("belldata #dim #Alice-keys #Bob-keys #trials\nNAlice and NBob > 3)\n");
        return 1;
        }
    key_choice NAlice = strtol(argv[2], NULL, 10);
    if(NAlice < 1)
        {
        printf("Enter 1 or more\n");
        return 1;
        }
    key_choice NBob = NAlice;
    if(argc > 3)
        NBob = strtol(argv[3], NULL, 10);
    if(NBob < 1)
        {
        printf("Enter 1 or more\n");
        return 1;
        }
    trial Ntrials = 0;
    if(argc > 4)
        Ntrials = strtol(argv[4], NULL, 10);
    printf("Ntrials %d\n", Ntrials);
    trial noff = 0;
    if(argc > 5)
        {
        noff = strtol(argv[5], NULL, 10);
        printf("Bob's counting is off by %d trials.\n", noff);
        }

    printf("Number of spatial dimensions %d\n", Dimensions);
    printf("Alice  has %d keys\n", NAlice);
    printf("Bob    has %d keys\n", NBob);
    printf("Number of trials: %d (Each trial consists of a pair of key choices and colour observations.)\n", Ntrials);
    printf("Total number of combinations of keys: %d\n", NAlice * NBob);

    doExperiment(Dimensions, NAlice, NBob, Ntrials, noff);
    return 0;
    }
