#include "../include/utilities.h"
#include "../include/ParametersManager.h"
#include "../include/ControlPickStrategy.h"
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>

extern "C" {

// ---------------------------
// Paramètres 
// ---------------------------
std::string paramsFile = "FarmAgroEco_params.json";

double rB = 0.4;
double KB = 8000.0;
double hB = 0.02;

double alphaS = 0.08;
double betaS = 0.0005;
double gammaS = 0.08;
double lambdaS = 0.05;

double pB = 0.5;           
double deltaC = 0.2;       
double rhoS = 0.01;
double cS = 20.0;
double cW = 15.0;

double R = 2.0;
double theta = 0.05;
double etaW = 0.02;
double deltaW = 0.01;

double kS = 100.0;
double kW = 100.0;
double kC = 100.0;
double Smax = 1.0;

// ---------------------------
// Dynamics 
// ---------------------------
void dynamics(const double *x, const double *u, double *image)
{
    double B = x[0];
    double S = x[1];
    double C = x[2];
    double W = x[3];

    double uB = u[0];
    double uS = u[1];
    double uW = u[2];

    double denomS = std::max(S + kS, 1e-8);
    double denomW = std::max(W + kW, 1e-8);
    double denomC = std::max(C + kC, 1e-8);

    double fS = S / denomS;
    double fW = W / denomW;

    double growthB = rB * B * (1.0 - B / KB);

    image[0] = growthB * fS * fW - (hB + uB) * B;

    image[1] = alphaS * (Smax - S)
               - betaS * B * S
               + gammaS * uS * (C / denomC)
               + lambdaS * fW;

    image[2] = pB * B
               - deltaC * C
               - rhoS * ((Smax - S) / Smax) * C
               - cS * uS
               - cW * uW;

    image[3] = R
               + uW
               + theta * uS * (C / denomC)
               - etaW * B * fW
               - deltaW * W;

    for (int i = 0; i < 4; ++i) {
        if (std::isnan(image[i]) || std::isinf(image[i])) image[i] = 0.0;
    }
}

// ---------------------------
// Contraintes d'état 
// ---------------------------
double constraintsX(const double *x)
{
    double B = x[0];
    double S = x[1];
    double C = x[2];
    double W = x[3];

    bool okB = (B >= 0.0 && B <= 8000.0);      
    bool okS = (S >= 0.0 && S <= 1.0);         
    bool okC = (C >= 0.0 && C <= 2000.0);      
    bool okW = (W >= 0.0 && W <= 500.0);       

    return (okB && okS && okC && okW) ? 1.0 : PLUS_INF;
}

// ---------------------------
// Contraintes état-commande 
// ---------------------------
double constraintsXU(const double *x, const double *u)
{
    if (constraintsX(x) == PLUS_INF) return PLUS_INF;

    bool okU0 = (u[0] >= -0.5 && u[0] <= 0.5);    
    bool okU1 = (u[1] >= -0.5 && u[1] <= 0.5);    
    bool okU2 = (u[2] >= -0.5 && u[2] <= 0.5);    

    return (okU0 && okU1 && okU2) ? 1.0 : PLUS_INF;
}

// ---------------------------
// Temporal control 
// ---------------------------
void temporalControl(double time, int, int, double *control)
{
    control[0] = 0.05;
    control[1] = 0.05;
    control[2] = 0.05;
}

// ---------------------------
// Jacobian analytique (4x4)
// ---------------------------
void jacobian(const double *x, const double *u, double ** jacob)
{
    double B = x[0];
    double S = x[1];
    double C = x[2];
    double W = x[3];

    double uB = u[0];
    double uS = u[1];
    double uW = u[2];

    double denomS = std::max(S + kS, 1e-8);
    double denomW = std::max(W + kW, 1e-8);
    double denomC = std::max(C + kC, 1e-8);

    double fS = S / denomS;
    double fW = W / denomW;

    double dfS_dS = kS / (denomS * denomS);
    double dfW_dW = kW / (denomW * denomW);
    double dCoverdC = kC / (denomC * denomC);

    double gB = rB * B * (1.0 - B / KB);
    double dgB_dB = rB * (1.0 - 2.0 * B / KB);

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            jacob[i][j] = 0.0;

    jacob[0][0] = dgB_dB * fS * fW - (hB + uB);
    jacob[0][1] = gB * dfS_dS * fW;
    jacob[0][2] = 0.0;
    jacob[0][3] = gB * fS * dfW_dW;

    jacob[1][0] = - betaS * S;
    jacob[1][1] = - alphaS - betaS * B;
    jacob[1][2] = gammaS * uS * dCoverdC;
    jacob[1][3] = lambdaS * dfW_dW;

    jacob[2][0] = pB;
    jacob[2][1] = rhoS * C / Smax;
    jacob[2][2] = - deltaC - rhoS * ((Smax - S) / Smax);
    jacob[2][3] = 0.0;

    jacob[3][0] = - etaW * fW;
    jacob[3][1] = 0.0;
    jacob[3][2] = theta * uS * dCoverdC;
    jacob[3][3] = - etaW * B * dfW_dW - deltaW;

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (std::isnan(jacob[i][j]) || std::isinf(jacob[i][j])) jacob[i][j] = 0.0;
}

// ---------------------------
// Bornes locales
// ---------------------------
void localDynBounds(const double *x, double * bound)
{
    bound[0] = 500.0;
    bound[1] = 0.5;
    bound[2] = 300.0;      
    bound[3] = 50.0;
}

} // extern "C"

// ---------------------------
// Target durcie
// ---------------------------
bool target(const double *x)
{
    double B = x[0];
    double S = x[1];
    double C = x[2];
    double W = x[3];

    bool okB = (B >= 500.0 && B <= 7500.0);    
    bool okS = (S >= 0.5 && S <= 0.95);      
    bool okC = (C >= 500.0 && C <= 1800.0);    
    bool okW = (W >= 50.0 && W <= 480.0);    

    return (okB && okS && okC && okW);
}

// ---------------------------
// Contrôle d'éligibilité pour trajectoire
// ---------------------------
extern "C" {
double controlEligibilityForTraj_fd(const double *x, const double *u) {
    return constraintsXU(x, u);
}
}

// ---------------------------
// Test rapide de target
// ---------------------------
int main() {
    double test1[4] = {600, 0.6, 600, 100};  // OK
    double test2[4] = {600, 0.4, 600, 100};  // S 
    double test3[4] = {600, 0.6, 400, 100};  // C 

    std::cout << "Test1: " << target(test1) << std::endl; // 1
    std::cout << "Test2: " << target(test2) << std::endl; // 0
    std::cout << "Test3: " << target(test3) << std::endl; // 0

    return 0;
}

