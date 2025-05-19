#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>


// Marouane
void libererMatrice(double** matrice, int lignes) {
    int i;
    for (i = 0; i < lignes; i++) {
        free(matrice[i]);
    }
    free(matrice);
}

void AfficherMatrice( double** matrice, int lignes, int colonnes) {
    int i,j;
    for (i = 0; i < lignes; i++) {
        printf("  [");
        for (j = 0; j < colonnes; j++) {
            printf("%10.6f ", matrice[i][j]);
        }
        printf("]\n");
    }
}

void AfficherVecteur( const double* vecteur, int taille) {
    int i;
    printf("  [ ");
    for ( i = 0; i < taille; i++) {
        printf("%10.6f ", vecteur[i]);
    }
    printf("]\n");
}

// Saad et Imaddine
double** CalculerMatriceRendements(double** matricePrix, int nombrePeriodes, int nombreActifs) {
    int i,j;
    double prixPrecedent, prixCourant;
    int nombrePeriodesRendement;

    nombrePeriodesRendement = nombrePeriodes - 1;
    
    double** matriceRendements = (double**)malloc((size_t)nombrePeriodesRendement * sizeof(double*));
    for (i = 0; i < nombrePeriodesRendement; i++) {
        matriceRendements[i] = (double*)malloc((size_t)nombreActifs * sizeof(double));
    }
    
    for (j = 0; j < nombreActifs; j++) {
        for (i = 0; i < nombrePeriodesRendement; i++) { // i correspond à la période de rendement
            prixPrecedent = matricePrix[i][j];
            prixCourant = matricePrix[i + 1][j];
            if (fabs(prixPrecedent) < 1e-9) {
                matriceRendements[i][j] = 0.0;
            } else {
                matriceRendements[i][j] = (prixCourant - prixPrecedent) / prixPrecedent;
            }
        }
    }
    return matriceRendements;
}

double* CalculerVecteurRendementsMoyens(double** matriceRendements, int nombrePeriodesRendement, int nombreActifs) {
    int i,j;
    double somme;
    double* vecteurRendementsMoyens = (double*)malloc((size_t)nombreActifs * sizeof(double));

    for (j = 0; j < nombreActifs; j++) {
        somme = 0.0;
        for (i = 0; i < nombrePeriodesRendement; i++) {
            somme += matriceRendements[i][j];
        }
        vecteurRendementsMoyens[j] = somme / (double)nombrePeriodesRendement;
    }
    return vecteurRendementsMoyens;
}

double** CalculerMatriceCovariance(double** matriceRendements, const double* vecteurRendementsMoyens, int nombrePeriodesRendement, int nombreActifs) {

    int i,j,k;
    double denominateur,somme,diff_i,diff_j;
    double** matriceCovariance = (double**)malloc((size_t)nombreActifs * sizeof(double*));
    for (i = 0; i < nombreActifs; i++) {
        matriceCovariance[i] = (double*)malloc((size_t)nombreActifs * sizeof(double));
    }

    denominateur = (double)(nombrePeriodesRendement - 1);

    for (i = 0; i < nombreActifs; i++) {
        for (j = 0; j <= i; j++) { // Calculer seulement la moitié inférieure (et diagonale)
            somme = 0.0;
            for (k = 0; k < nombrePeriodesRendement; k++) {
                diff_i = matriceRendements[k][i] - vecteurRendementsMoyens[i];
                diff_j = matriceRendements[k][j] - vecteurRendementsMoyens[j];
                somme += diff_i * diff_j;
            }
            matriceCovariance[i][j] = somme / denominateur;
            if (i != j) { // Remplir la partie supérieure par symétrie
                matriceCovariance[j][i] = matriceCovariance[i][j];
            }
        }
    }
    return matriceCovariance;
}

double** ConstruireMatriceKKT(double** matriceCovariance, const double* vecteurRendementsMoyens, int nombreActifs) {

    int i,j;
    int tailleKKT = nombreActifs + 2;

    double** kkt = (double**)malloc((size_t)tailleKKT * sizeof(double*));
    for (int i = 0; i < tailleKKT; i++) {
        kkt[i] = (double*)calloc((size_t)tailleKKT, sizeof(double)); // calloc initialise à zéro
    }

    // Bloc 1: Matrice de covariance
    for ( i = 0; i < nombreActifs; i++) {
        for (int j = 0; j < nombreActifs; j++) {
            kkt[i][j] = matriceCovariance[i][j];
        }
    }

    // Bloc 2 et 3: Vecteur des rendements moyens et vecteur de uns
    for ( i = 0; i < nombreActifs; i++) {
        kkt[i][nombreActifs]     = vecteurRendementsMoyens[i];
        kkt[i][nombreActifs + 1] = 1.0;

        kkt[nombreActifs][i]     = vecteurRendementsMoyens[i];
        kkt[nombreActifs + 1][i] = 1.0;
    }

    // Le bloc inférieur droit (2x2) reste à zéro grâce à calloc :)
    return kkt;
}

// Abderrahmane et Faycal
double* ResoudreSystemeLineaire(double** matriceA_orig, const double* Y_orig, int N_eq) {  

    double** LU = (double**)malloc((size_t)N_eq * sizeof(double*));
    int* P = (int*)malloc((size_t)N_eq * sizeof(int));
    int i,j,k,i_p,tempP;
    double maxValAbs,somme;
    double* temp_col_ptr;
    double* X = (double*)malloc((size_t)N_eq * sizeof(double));
    double* Y_P = (double*)malloc((size_t)N_eq * sizeof(double));

    // Vérification de l'allocation mémoire et copie de la matrice A en LU-------
    if (matriceA_orig == NULL || Y_orig == NULL) {
        return NULL;
    }

    if (LU == NULL || X == NULL || Y_P == NULL || P == NULL) {
        if (LU) { for (i = 0; i < N_eq; ++i) free(LU[i]); free(LU); } 
        free(X); free(Y_P); free(P);
        return NULL;
    }

    for (i = 0; i < N_eq; i++) {
        LU[i] = (double*)malloc((size_t)N_eq * sizeof(double));
        if (LU[i] == NULL) {
            libererMatrice(LU, i);
            free(X); free(Y_P); free(P);
            return NULL;
        }
        if (matriceA_orig[i] == NULL) {
            libererMatrice(LU, i + 1);
            free(X); free(Y_P); free(P);
            return NULL;
        }
        memcpy(LU[i], matriceA_orig[i], (size_t)N_eq * sizeof(double)); // copie une ligne entière de la matrice matriceA_orig dans la matrice LU
    }
    // ---------------------------------------------------------------------------

    // Initialisation du vecteur de permutation
    for (i = 0; i < N_eq; i++) P[i] = i;

    // Décomposition LU avec pivotage partiel
    for (k = 0; k < N_eq; k++) {
        i_p = k; // Ligne du pivot
        maxValAbs = fabs(LU[k][k]);

        for (i = k + 1; i < N_eq; i++) {
            if (fabs(LU[i][k]) > maxValAbs) {
                maxValAbs = fabs(LU[i][k]);
                i_p = i;
            }
        }

        // Permutation des lignes si nécessaire
        if (i_p != k) {
            temp_col_ptr = LU[k];
            LU[k] = LU[i_p];
            LU[i_p] = temp_col_ptr;

            tempP = P[k];
            P[k] = P[i_p];
            P[i_p] = tempP;
        }

        // Vérification du pivot
        if (fabs(LU[k][k]) < 1e-12) {
            fprintf(stderr, "    ResoudreSystemeLineaire: ATTENTION - Pivot proche de zéro à LU[%d][%d] = %e. Matrice singulière ou mal conditionnée.\n", k, k, LU[k][k]);
        }

        // Calcul des éléments de L et U
        for (i = k + 1; i < N_eq; i++) {
            if (fabs(LU[k][k]) > 1e-12) {
                 LU[i][k] = LU[i][k] / LU[k][k]; // Éléments de L
            } else {
                 LU[i][k] = 0.0;
            }
            for (j = k + 1; j < N_eq; j++) {
                LU[i][j] -= LU[i][k] * LU[k][j]; // Éléments de U
            }
        }
    }

    // Résolution de Ly = Pb (Substitution avant)
    for (i = 0; i < N_eq; i++) Y_P[i] = Y_orig[P[i]];

    for (i = 0; i < N_eq; i++) {
        somme = Y_P[i];
        for (j = 0; j < i; j++) {
            somme -= LU[i][j] * Y_P[j];
        }
        Y_P[i] = somme;
    }

    // Résolution de Ux = y (Substitution arrière)
    for (i = N_eq - 1; i >= 0; i--) {
        somme = Y_P[i];
        for (j = i + 1; j < N_eq; j++) {
            somme -= LU[i][j] * X[j];
        }
        if (fabs(LU[i][i]) < 1e-12) {
            X[i] = 0.0;
        } else {
            X[i] = somme / LU[i][i];
        }
    }

    libererMatrice(LU, N_eq);
    free(P);
    free(Y_P);
    return X;
}

double* ResoudreSystemeOptimisation(double** matriceKKT, double rendementCible, int nombreActifs) {

    int tailleKKT = nombreActifs + 2;
    double* Y = (double*)calloc((size_t)tailleKKT, sizeof(double));
    double* poids_finaux = (double*)malloc((size_t)nombreActifs * sizeof(double));
    double* X_solution_iteration = NULL;

    const int MAX_ITERATIONS_OPTIM = nombreActifs; // Pour que la boucle s'exécute une fois comme l'original
    int iteration_optim_count;
    bool contrainte_negative_activee;

    Y[nombreActifs]     = rendementCible;
    Y[nombreActifs + 1] = 1.0;

    iteration_optim_count = 0;
    do {
        iteration_optim_count++;
        free(X_solution_iteration); // Libérer la solution de l'itération précédente s'il y en a une

        // Résoudre le système KKT*X = Y pour obtenir les poids bruts et les multiplicateurs de Lagrange
        X_solution_iteration = ResoudreSystemeLineaire(matriceKKT, Y, tailleKKT);

        contrainte_negative_activee = false;

        // Troncature des poids négatifs (simplification, pas une vraie méthode d'ensemble actif)
        for (int i = 0; i < nombreActifs; i++) {
            if (X_solution_iteration[i] < -1e-9) { // Utiliser une petite tolérance
                X_solution_iteration[i] = 0.0;
                contrainte_negative_activee = true; // Non utilisé car MAX_ITERATIONS_OPTIM = 1
            }
        }
    } while (contrainte_negative_activee && iteration_optim_count < MAX_ITERATIONS_OPTIM);

    // Normalisation des poids pour qu'ils somment à 1
    double somme_des_poids = 0.0;
    for (int i = 0; i < nombreActifs; i++) {
        somme_des_poids += X_solution_iteration[i];
    }

    if (fabs(somme_des_poids) < 1e-9) {
        for (int i = 0; i < nombreActifs; i++){
            poids_finaux[i] = 1.0 / (double)nombreActifs;
        }
    } else {
        for (int i = 0; i < nombreActifs; i++){
            poids_finaux[i] = X_solution_iteration[i] / somme_des_poids;
        }
    }

    free(Y);
    free(X_solution_iteration);
    return poids_finaux;
}

double CalculerRendementPortefeuille(const double* poids, const double* rendementsMoyens, int nombreActifs) {

    double rendementPortefeuille = 0.0;
    for (int i = 0; i < nombreActifs; i++) {
        rendementPortefeuille += poids[i] * rendementsMoyens[i];
    }
    return rendementPortefeuille;
}

double CalculerVolatilitePortefeuille(const double* poids, double** matriceCovariance, int nombreActifs) {

    double variancePortefeuille = 0.0;
    for (int i = 0; i < nombreActifs; i++) {
        for (int j = 0; j < nombreActifs; j++) {
            variancePortefeuille += poids[i] * poids[j] * matriceCovariance[i][j];
        }
    }
    if (variancePortefeuille < 0) {
        variancePortefeuille = 0.0;
    }
    return sqrt(variancePortefeuille);
}

int main() {

    printf("Debut du programme d'optimisation de portefeuille.\n");

    const int nombreJours = 89;
    const int nombreActifs = 7;
    const double rendementCible = 0.0005;
    
    printf("Initialisation des donnees: nombreJours=%d, nombreActifs=%d\n", nombreJours, nombreActifs);

    // Allocation pour matricePrix
    double** matricePrix = (double**)malloc((size_t)nombreJours * sizeof(double*));
    if (matricePrix == NULL) { return 1; }
    for (int i = 0; i < nombreJours; i++) {
        matricePrix[i] = (double*)malloc((size_t)nombreActifs * sizeof(double));
        if (matricePrix[i] == NULL) {
            libererMatrice(matricePrix, i);
            return 1;
        }
    }


    const char* nomsActifs[] = {"NVDA", "TSLA", "BA", "AAPL", "BABA", "GOOGL", "AMZN"};
    // Données de prix (copiées depuis votre exemple)
    double tempPrix[89][7] = {
        //NVDA    TSLA     BA     AAPL   BABA    GOOGL   AMZN
        {138.31, 379.28, 171.87, 243.85, 84.95, 189.43, 220.22},
        {144.47, 410.44, 169.90, 243.36, 85.54, 191.79, 224.19},
        {149.43, 411.05, 170.78, 245.00, 85.52, 196.87, 227.61},
        {140.14, 394.36, 172.51, 242.21, 84.48, 195.49, 222.11},
        {140.11, 394.94, 171.76, 242.70, 83.69, 193.95, 222.13},
        {135.91, 394.74, 172.00, 236.85, 80.53, 192.04, 218.94},
        {133.23, 403.31, 170.57, 234.40, 80.54, 191.01, 218.46},
        {131.76, 396.36, 167.02, 233.28, 81.68, 189.66, 217.76},
        {136.24, 428.22, 166.20, 237.87, 82.44, 195.55, 223.35},
        {133.57, 413.82, 168.93, 228.26, 82.43, 192.91, 220.66},
        {137.71, 426.50, 171.09, 229.98, 85.12, 196.00, 225.94},
        {140.83, 424.07, 175.56, 222.64, 85.38, 198.05, 230.71},
        {147.07, 415.11, 174.80, 223.83, 86.40, 198.37, 235.01},
        {147.22, 412.38, 178.50, 223.66, 86.10, 197.98, 235.42},
        {142.62, 406.58, 176.06, 222.78, 89.14, 200.21, 234.85},
        {118.42, 397.15, 175.16, 229.86, 89.99, 191.81, 235.42},
        {128.99, 398.09, 177.78, 238.26, 96.03, 195.30, 238.15},
        {123.70, 389.10, 173.66, 239.36, 96.72, 195.41, 237.07},
        {124.65, 400.28, 179.53, 237.59, 102.74, 200.87, 234.64},
        {120.07, 404.60, 176.52, 236.00, 98.84, 204.02, 237.68},
        {116.66, 383.68, 175.87, 228.01, 98.61, 201.23, 237.42},
        {118.65, 392.21, 176.23, 232.80, 102.35, 206.38, 242.06},
        {124.83, 378.17, 181.84, 232.47, 99.28, 191.33, 236.17},
        {128.68, 374.32, 184.80, 233.22, 100.38, 191.60, 238.83},
        {129.84, 361.62, 181.49, 227.63, 103.51, 185.34, 229.15},
        {133.57, 350.73, 180.55, 227.65, 111.32, 186.47, 233.14},
        {132.80, 328.50, 180.44, 232.62, 112.78, 185.32, 232.76},
        {131.14, 336.51, 186.25, 236.87, 118.33, 183.61, 228.93},
        {135.29, 355.94, 185.44, 241.53, 119.54, 186.14, 230.37},
        {138.85, 355.84, 184.42, 244.60, 124.73, 185.23, 228.68},
        {139.40, 354.11, 184.97, 244.47, 126.90, 183.77, 226.65},
        {139.23, 360.56, 186.15, 244.87, 125.79, 185.27, 226.63},
        {140.11, 354.40, 180.88, 245.83, 135.97, 184.56, 222.88},
        {134.43, 337.80, 177.15, 245.55, 143.75, 179.66, 216.58},
        {130.28, 330.53, 179.91, 247.10, 129.04, 179.25, 212.71},
        {126.63, 302.80, 178.27, 247.04, 134.01, 175.42, 212.80},
        {131.28, 290.80, 173.04, 240.36, 139.08, 172.73, 214.35},
        {120.15, 281.95, 173.83, 237.30, 136.55, 168.50, 208.74},
        {124.92, 292.98, 174.63, 241.84, 132.51, 170.28, 212.28},
        {114.06, 284.65, 170.06, 238.03, 130.81, 167.01, 205.02},
        {115.99, 272.04, 158.90, 235.93, 129.85, 170.92, 203.80},
        {117.30, 279.10, 163.16, 235.74, 141.03, 173.02, 208.36},
        {110.57, 263.45, 158.43, 235.33, 139.95, 172.35, 200.70},
        {112.69, 262.67, 154.18, 239.07, 140.62, 173.86, 199.25},
        {106.98, 222.15, 148.15, 227.48, 132.54, 165.87, 194.54},
        {108.76, 230.58, 154.06, 220.84, 139.02, 164.04, 196.59},
        {115.74, 248.09, 158.80, 216.98, 137.14, 167.11, 198.89},
        {115.58, 240.68, 159.32, 209.68, 138.35, 162.76, 193.89},
        {121.67, 249.98, 161.81, 213.49, 141.10, 165.49, 197.95},
        {119.53, 238.01, 161.85, 214.00, 147.57, 164.29, 195.74},
        {115.43, 225.31, 161.57, 212.69, 142.74, 160.67, 192.82},
        {117.52, 235.86, 172.62, 215.24, 143.20, 163.89, 195.54},
        {118.53, 236.26, 172.83, 214.10, 136.91, 162.80, 194.95},
        {117.70, 248.71, 178.11, 218.27, 135.14, 163.99, 196.21},
        {121.41, 278.39, 180.90, 220.73, 134.48, 167.68, 203.26},
        {120.69, 288.14, 182.59, 223.75, 132.75, 170.56, 205.71},
        {113.76, 272.06, 178.55, 221.53, 132.24, 165.06, 201.13},
        {111.43, 273.13, 179.11, 223.85, 135.63, 162.24, 201.36},
        {109.67, 263.55, 173.31, 217.90, 132.43, 154.33, 192.72},
        {108.38, 259.16, 170.55, 222.13, 132.23, 154.64, 190.26},
        {110.15, 268.46, 168.17, 223.19, 132.70, 157.07, 192.17},
        {110.42, 282.76, 168.56, 223.89, 129.79, 157.04, 196.01},
        {101.80, 267.28, 150.91, 203.19, 129.33, 150.72, 178.41},
        {94.31, 239.43, 136.59, 188.38, 116.54, 145.60, 171.00},
        {97.64, 233.29, 138.86, 181.46, 105.98, 146.75, 175.26},
        {96.30, 221.86, 139.39, 172.42, 99.37, 144.70, 170.66},
        {114.33, 272.20, 160.82, 198.85, 104.78, 158.71, 191.10},
        {107.57, 252.40, 155.52, 190.42, 104.18, 152.82, 181.22},
        {110.93, 252.31, 156.84, 198.15, 107.73, 157.14, 184.87},
        {110.71, 252.35, 159.28, 202.52, 113.97, 159.07, 182.12},
        {112.20, 254.11, 155.52, 202.14, 112.28, 156.31, 179.59},
        {104.49, 241.55, 156.47, 194.27, 106.75, 153.33, 174.33},
        {101.49, 241.37, 161.90, 196.98, 108.87, 151.16, 172.61},
        {96.91, 227.50, 159.34, 193.16, 110.15, 147.67, 167.32},
        {98.89, 237.97, 162.52, 199.74, 115.88, 151.47, 173.18},
        {102.71, 250.74, 172.37, 204.60, 118.97, 155.35, 180.60},
        {106.43, 259.51, 176.26, 208.37, 119.29, 159.28, 186.54},
        {111.01, 284.95, 177.95, 209.28, 120.28, 161.96, 188.99},
        {108.73, 285.88, 182.30, 210.14, 118.37, 160.61, 187.70},
        {109.02, 292.03, 182.00, 211.21, 118.88, 160.16, 187.39},
        {108.92, 282.16, 183.24, 212.50, 119.43, 158.80, 184.42},
        {111.61, 280.52, 182.89, 213.32, 120.53, 161.30, 190.20},
        {114.50, 287.21, 185.46, 205.35, 125.76, 164.03, 189.98},
        {113.82, 280.26, 186.46, 198.89, 126.57, 164.21, 186.35},
        {113.54, 275.35, 185.96, 198.51, 127.66, 163.23, 185.01},
        {117.06, 276.22, 185.56, 196.25, 123.23, 151.38, 188.71},
        {117.37, 284.82, 191.70, 197.49, 125.79, 154.28, 192.08},
        {116.65, 298.26, 194.85, 198.53, 125.33, 152.75, 193.06}
    };

    for(int i=0; i<nombreJours; ++i) {
        for(int j=0; j<nombreActifs; ++j) {
            matricePrix[i][j] = tempPrix[i][j];
        }
    }

    printf("\n--- ETAPE 1: Calcul Matrice Rendements ---\n");
    double** matriceRendements = CalculerMatriceRendements(matricePrix, nombreJours, nombreActifs);
    //AfficherMatrice( matriceRendements, nombreJours - 1, nombreActifs);

    printf("\n--- ETAPE 2: Calcul Vecteur Rendements Moyens ---\n");
    int nombrePeriodesRendement = nombreJours - 1;
    double* vecteurRendementsMoyens = CalculerVecteurRendementsMoyens(matriceRendements, nombrePeriodesRendement, nombreActifs);
    AfficherVecteur( vecteurRendementsMoyens, nombreActifs);

    printf("\n--- ETAPE 3: Calcul Matrice Covariance ---\n");
    double** matriceCovariance = NULL;
    if (nombrePeriodesRendement >= 1) {
        matriceCovariance = CalculerMatriceCovariance(matriceRendements, vecteurRendementsMoyens, nombrePeriodesRendement, nombreActifs);
        AfficherMatrice( matriceCovariance, nombreActifs, nombreActifs);
    }

    printf("\n--- ETAPE 4: Construction Matrice KKT ---\n");
    double** matriceKKT = ConstruireMatriceKKT(matriceCovariance, vecteurRendementsMoyens, nombreActifs);
    AfficherMatrice( matriceKKT, nombreActifs + 2, nombreActifs + 2);

    printf("\n--- ETAPE 5: Resolution Systeme Optimisation ---\n");
    printf("main: Rendement cible par periode = %.6f (%.4f %%)\n", rendementCible, rendementCible* 100.0);
    double* poidsOptimaux = ResoudreSystemeOptimisation(matriceKKT, rendementCible, nombreActifs);

    printf("\n--- RESULTATS FINALS ---\n");
    printf("Poids optimaux (format liste) :\n");
    for (int i = 0; i < nombreActifs; i++) {
        printf("  Actif %-5s : %10.6f (%8.4f %%)\n", nomsActifs[i], poidsOptimaux[i], poidsOptimaux[i]*100.0);
    }

    // Calcul du rendement avant optimisation avec des poids égaux
    double* poidsEgaux = (double*)malloc((size_t)nombreActifs * sizeof(double));

    for (int i = 0; i < nombreActifs; i++) {
        poidsEgaux[i] = 1.0 / nombreActifs;

        
    }

    // Calcul de la différence relative entre le rendement avant optimisation

    double rendementPtfAvant = CalculerRendementPortefeuille(poidsEgaux, vecteurRendementsMoyens, nombreActifs);
    double volatilitePtfAvant = CalculerVolatilitePortefeuille(poidsEgaux, matriceCovariance, nombreActifs);

    printf("\nRendement avant optimisation du portefeuille [Chaque poids est (1/Nombre actifs)]: %10.6f (%8.4f %% par periode)\n", rendementPtfAvant, rendementPtfAvant * 100.0);
    printf("Volatilite avant optimisation du portefeuille [Chaque poids est (1/Nombre actifs)]: %10.6f (%8.4f %% par periode)\n", volatilitePtfAvant, volatilitePtfAvant * 100.0);

    double rendementPtfOptimise = CalculerRendementPortefeuille(poidsOptimaux, vecteurRendementsMoyens, nombreActifs);
    double volatilitePtfOptimise= CalculerVolatilitePortefeuille(poidsOptimaux, matriceCovariance, nombreActifs);

    printf("\nRendement apres optimisation du portefeuille [les poids sont optimaux]: %10.6f (%8.4f %% par periode)\n", rendementPtfOptimise, rendementPtfOptimise * 100.0);
    printf("Volatilite apres optimisation du portefeuille [les poids sont optimaux]: %10.6f (%8.4f %% par periode)\n", volatilitePtfOptimise, volatilitePtfOptimise * 100.0);

    
    double differenceRelativeRendement = (rendementPtfOptimise - rendementPtfAvant) / fabs(rendementPtfAvant) * 100.0;

    // Calcul de la différence relative entre la volatilité avant et après optimisation
    double differenceRelativeVolatilite = (volatilitePtfAvant - volatilitePtfOptimise) / fabs(volatilitePtfAvant) * 100.0;

    printf("\nAmelioration relative du rendement apres optimisation : %.4f %%\n", differenceRelativeRendement);
    printf("Amelioration relative de la volatilite apres optimisation : %.4f %%\n", differenceRelativeVolatilite);

    printf("\nLiberation de la memoire...\n");
    libererMatrice(matricePrix, nombreJours);
    libererMatrice(matriceRendements, nombrePeriodesRendement);
    libererMatrice(matriceCovariance, nombreActifs);
    libererMatrice(matriceKKT, nombreActifs + 2);
    free(vecteurRendementsMoyens);
    free(poidsOptimaux);
    free(poidsEgaux);

    printf("Fin du programme.\n");
    return 0;
}