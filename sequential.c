#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This functions finds the determinant of Matrix
double determinantOfMatrix(double mat[3][3])
{
    double ans;
    ans = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2])
          - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
          + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    return ans;
}

// This function finds the solution of system of
// linear equations using cramer's rule
int findSolution(double coeff[3][4])
{
    double x, y, z;
    // Matrix d using coeff as given in cramer's rule
    double d[3][3] = {
        { coeff[0][0], coeff[0][1], coeff[0][2] },
        { coeff[1][0], coeff[1][1], coeff[1][2] },
        { coeff[2][0], coeff[2][1], coeff[2][2] },
    };
    // Matrix d1 using coeff as given in cramer's rule
    double d1[3][3] = {
        { coeff[0][3], coeff[0][1], coeff[0][2] },
        { coeff[1][3], coeff[1][1], coeff[1][2] },
        { coeff[2][3], coeff[2][1], coeff[2][2] },
    };
    // Matrix d2 using coeff as given in cramer's rule
    double d2[3][3] = {
        { coeff[0][0], coeff[0][3], coeff[0][2] },
        { coeff[1][0], coeff[1][3], coeff[1][2] },
        { coeff[2][0], coeff[2][3], coeff[2][2] },
    };
    // Matrix d3 using coeff as given in cramer's rule
    double d3[3][3] = {
        { coeff[0][0], coeff[0][1], coeff[0][3] },
        { coeff[1][0], coeff[1][1], coeff[1][3] },
        { coeff[2][0], coeff[2][1], coeff[2][3] },
    };

    // Calculating Determinant of Matrices d, d1, d2, d3
    double D = determinantOfMatrix(d);
    double D1 = determinantOfMatrix(d1);
    double D2 = determinantOfMatrix(d2);
    double D3 = determinantOfMatrix(d3);

    // Case 1
    if (D != 0) {
        // Coeff have a unique solution. Apply Cramer's Rule
        x = D1 / D;
        y = D2 / D;
        z = D3 / D; // calculating z using cramer's rule

        return x;
    }
    // Case 2
    else {
        if (D1 == 0 && D2 == 0 && D3 == 0)
            printf("Infinite solutions\n");
        else if (D1 != 0 || D2 != 0 || D3 != 0)
            printf("No solutions\n");
    }
}

int main(){
    // Define inputs
    double W, Height, Length, r, Wr, StrainRate, dH;
    double H1, H2, H3;
    W = 8.7;

    // Initial dimensions of model
    Height = 0.05;
    Length = 0.09;

    // Roller properties
    r = 0.0026;
    Wr = 0.035;

    StrainRate = 0.00125;

    // Considering one minute of time has passed
    dH = StrainRate;
    double theta = atan(dH/Height);

    // Iterating over top layer
    printf("---------------------Start of Top layer iteration---------------------\n");
    int j;
    H3 = 0;
    for(j = 1; j <= 17; j++){
        double top[3][4] = {
            {0, -(1+cos(theta)), 0, (Wr-W)},
            {1, sin(theta), -2, 0},
            {17*r, -(((2*j-1)*r+16*r*sin(theta))+cos(theta)*((2*j-1)*r+15*r*sin(theta))), -34*r, (Wr-W)*(2*j-1)*r+16*r*sin(theta)}
        };

        H1 = H1 + findSolution(top);
        printf("\n");
    }
    printf("---------------------End of Top layer iteration---------------------\n");

    // Iterating over intermediate layers
    printf("---------------------Start of Intermediate layer iteration---------------------\n");
    int i;
    double Wij;
    for(i = 1; i <= 6; i++){
        for(j = 1; j <= 17; j++){
            Wij = W + (8-i)*Wr;
            double inter[3][4] = {
                {0, -2*cos(theta), 0, -Wij*sin(theta)},
                {1, 2*sin(theta), -2, Wij*sin(theta)},
                {(2*i+1)*r, r*sin(theta)*(4*i + 1 + cos(theta))-2*cos(theta)*((2*j-1)*r+2*i*r*sin(theta)), -2*(2*i+1)*r, Wij*(2*(i+1))*r*sin(theta)-cos(theta)*((2*j-1)*r+(2*i+1)*r*sin(theta)) - Wr*((2*j-1)*r+2*i*r*sin(theta))}
            };

            H2 = H2 + findSolution(inter);
            printf("\n");
        }
    }
    printf("---------------------End of Intermediate layer iteration---------------------\n");

    // Iterating over bottom layer
    printf("---------------------Start of Bottom layer iteration---------------------\n");
    double Wj;
    for (j = 1; j <= 17; j++) {
        Wj = W + 8*Wr;
        double bottom[3][4] = {
            {0. -(1+cos(theta)), 0, -(Wr + Wj*cos(theta))},
            {1, sin(theta), -2, Wj},
            {r, -r*((2*j-1)*(1+cos(theta))+sin(theta)), -2*r, Wj*r*(sin(theta) - cos(theta)*((2*j-1)+sin(theta))) - W*r}
        };

        H3 = H3 + findSolution(bottom);
        printf("\n");
    }
    printf("---------------------End of Bottom layer iteration---------------------\n");

    // Calculating total Horizontal force F
    double F = H1+H2+H3;
    printf("Total Horizontal force F = %f\n", F);
}
