#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

void multiAx(int m, int n, const float A[m][n],const float x[n],float y[m]){
    for (int i = 0; i < m; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            sum += A[i][j] * x[j];
        y[i] = (float)sum;
    }
}

void multiATx(int m, int n,
                       const float A[m][n],
                       const float x[m],
                       float y[n])
{
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int i = 0; i < m; i++)
            sum += A[i][j] * x[i];
        y[j] = (float)sum;
    }
}

 float norm(const float *v, int len) {
    double s = 0.0;
    for (int i = 0; i < len; i++) s += (double)v[i]*v[i];
    return (float)sqrt(s);
}

 void unit(float *v, int p) {
    float n = norm(v, p);
    if (n > 0.0f) {
        float inv = 1.0f / n;
        for (int i = 0; i < p; i++) v[i] *= inv;
    }
}

unsigned char alter(float x) {
    if (x < 0.0f)   return 0;
    if (x > 255.0f) return 255;
    return (unsigned char)(x + 0.5f);
}
 float body(int m, int n,const float A[m][n],
                                  float *u, float *v,
                                  float *temp1, float *temp2,
                                  int max, float tol)
{
    for (int j = 0; j < n; ++j) v[j] = (float)(j+2);
    unit(v, n);
    float oldsigma = 0.0f;
    for (int i = 0; i < max; ++i) {
        multiAx(m, n, A, v, temp1);
        float sigma = norm(temp1, m);
        if (sigma == 0.0f) return 0.0f;
        for (int i = 0; i < m; ++i) u[i] = temp1[i] / sigma;
        multiATx(m, n, A, u, temp2);
        unit(temp2, n);
        for (int j = 0; j < n; ++j) v[j] = temp2[j];
        multiAx(m, n, A, v, temp1);
        float newsigma = norm(temp1, m);
        if (fabsf(newsigma- oldsigma) < tol * newsigma) break;
        oldsigma =newsigma;
    }
    multiAx(m, n, A, v, temp1);
    float sigma = norm(temp1, m);
    if (sigma > 0.0f)
        for (int i = 0; i < m; ++i) u[i] = temp1[i] / sigma;
    return sigma;
}


 void add(int m, int n, float R[m][n],
                      const float *u, const float *v, float s)
{
    for (int i = 0; i < m; ++i) {
        float ui = u[i] * s;
        for (int j = 0; j < n; ++j)
            R[i][j] += ui * v[j];
    }
}

 void sub(int m, int n, float A[m][n],
                    const float *u, const float *v, float s)
{
    for (int i = 0; i < m; ++i) {
        float ui = u[i] * s;
        for (int j = 0; j < n; ++j)
            A[i][j] -= ui * v[j];
    }
}

double error(int m, int n,const float A[m][n], const float B[m][n]) {
    double sum = 0.0;
    double normA=0.0;
for(int j=0;j<m;j++)
    for(int i = 0; i < n; i++) {
        normA += (double)A[j][i]*(double)A[j][i];
        double diff = (double)A[j][i] - (double)B[j][i];
        sum += diff * diff;
    }
         double  a=(sqrt(sum))/sqrt(normA);
    return  a*100;
}

int main() {

    char input[200];
    printf("Enter input image filename: ");
    scanf("%199s", input);

    int k;
    printf("Enter k: ");
    scanf("%d", &k);

    char output[200];
     printf("Enter output image filename: ");
    scanf("%199s", output);


    int w, h, ch;
    unsigned char *img = stbi_load(input, &w, &h, &ch, 1);
    if (!img) return 1;

    int m = h, n = w;

    float (*A)[n]  = malloc(m*n*sizeof(float));
    float (*Aw)[n] = malloc(m*n*sizeof(float));
    float (*R)[n]  = calloc(m, n*sizeof(float));

    float *u  = malloc(m*sizeof(float));
    float *v  = malloc(n*sizeof(float));
    float *t1 = malloc(m*sizeof(float));
    float *t2 = malloc(n*sizeof(float));

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            A[i][j] = img[i*n + j];

    memcpy(Aw, A, m*n*sizeof(float));

    for (int t = 0; t < k; ++t) {
        float s = body(m, n,Aw, u, v, t1, t2, 200, 1e-4f);
        if (s < 1e-6f) break;
        add(m, n, R, u, v, s);
        sub(m, n, Aw, u, v, s);
    }
 
   double err=error(m,n,A,R);
   printf("Relative Error = %.2lf%%\n",err);

    unsigned char *out = malloc(m*n);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            out[i*n + j] = alter(R[i][j]);

    stbi_write_jpg(output, n, m, 1, out, 90);

    stbi_image_free(img);
    free(A); free(Aw); free(R);
    free(u); free(v); free(t1); free(t2);
    free(out);

    return 0;
}
