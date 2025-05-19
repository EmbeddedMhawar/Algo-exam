# Optimiseur Quadratique Contraint en C

Ce projet C implémente un solveur pour des problèmes d'optimisation quadratique sous contraintes d'égalité linéaires. Il est initialement conçu pour l'optimisation de portefeuille (trouver des poids d'actifs optimaux), mais le solveur principal peut être adapté à d'autres domaines.

## Description du Code

Le programme effectue les étapes suivantes :

1.  **Préparation des Données (Spécifique au Portefeuille) :**
    *   Charge ou génère une matrice de prix historiques d'actifs.
    *   Calcule la matrice des rendements périodiques à partir des prix.
    *   Calcule le vecteur des rendements moyens de chaque actif.
    *   Calcule la matrice de covariance des rendements des actifs.

2.  **Construction du Système d'Optimisation :**
    *   Formule le problème d'optimisation en utilisant les conditions de Karush-Kuhn-Tucker (KKT).
    *   Construit la matrice KKT qui intègre la fonction objectif (variance du portefeuille à minimiser), les contraintes d'égalité (rendement cible du portefeuille et somme des poids égale à 1).
    *   Construit le vecteur de droite (RHS) du système KKT.

3.  **Résolution du Système d'Optimisation :**
    *   Le cœur du programme est la fonction `ResoudreSystemeOptimisation`.
    *   Cette fonction prend la matrice KKT, un rendement cible et le nombre d'actifs.
    *   Elle appelle `ResoudreSystemeLineaire` pour résoudre le système `KKT * X = Y`, où `X` contient les poids des actifs et les multiplicateurs de Lagrange.
    *   **Gestion des Contraintes d'Inégalité (Simplifiée) :** La version actuelle gère la contrainte de non-négativité des poids (`w_i >= 0`) de manière simplifiée : après avoir résolu le système KKT (qui ne prend pas directement en compte cette inégalité), elle tronque tout poids négatif à zéro.
    *   Les poids résultants sont ensuite normalisés pour s'assurer que leur somme est exactement égale à 1.

4.  **Calcul des Métriques du Portefeuille :**
    *   Calcule le rendement attendu et la volatilité (écart-type) du portefeuille optimisé.
    *   Affiche les résultats.

## Fonctionnement du Solveur d'Optimisation

Le solveur d'optimisation au cœur de ce code est basé sur la résolution d'un système d'équations linéaires dérivé des conditions KKT.

### Problème d'Optimisation de Portefeuille (Exemple)

L'objectif typique est de minimiser la variance du portefeuille (risque) `σ_p^2` sous certaines contraintes :
*   **Minimiser :** `w^T * Σ * w`
    *   `w` : vecteur des poids des actifs
    *   `Σ` : matrice de covariance des rendements des actifs
*   **Sous les contraintes :**
    1.  `E^T * w = μ_cible` (Le rendement attendu du portefeuille est égal à un rendement cible `μ_cible`)
        *   `E` : vecteur des rendements moyens attendus des actifs
    2.  `1^T * w = 1` (La somme des poids des actifs est égale à 1)
    3.  `w_i >= 0` (Pas de vente à découvert – les poids sont non-négatifs). *Cette contrainte est gérée par troncature dans la version actuelle.*

### Conditions KKT et Système Linéaire

Pour les contraintes d'égalité, on introduit des multiplicateurs de Lagrange (λ). Le Lagrangien est :
`L(w, λ_1, λ_2) = w^T * Σ * w - λ_1 * (E^T * w - μ_cible) - λ_2 * (1^T * w - 1)`

En prenant les dérivées partielles par rapport à `w`, `λ_1`, et `λ_2` et en les égalant à zéro, on obtient un système d'équations linéaires :
`[ 2Σ   E   1 ] [  w  ]   [   0    ]`
`[ E^T  0   0 ] [ λ_1 ] = [ μ_cible ]`
`[ 1^T  0   0 ] [ λ_2 ]   [   1    ]`

*   La matrice `2Σ` peut être remplacée par `Σ` si la fonction objectif est formulée comme `(1/2)w^T * Σ * w`. Le code actuel utilise `Σ` dans `ConstruireMatriceKKT`.
*   Le vecteur `[0 ... 0]^T` dans le RHS a la taille du nombre d'actifs.

La fonction `ConstruireMatriceKKT` assemble la matrice de gauche, et `ResoudreSystemeOptimisation` prépare le vecteur de droite `Y_rhs`.

### `ResoudreSystemeLineaire`

Cette fonction résout un système d'équations linéaires `A*x = b` (ici, `KKT * X = Y_rhs`). L'implémentation actuelle utilise une **décomposition LU avec pivotage partiel**.
1.  **Décomposition LU :** La matrice `A` (ici `KKT`) est décomposée en un produit de deux matrices : une matrice triangulaire inférieure `L` et une matrice triangulaire supérieure `U`, de sorte que `P*A = L*U`, où `P` est une matrice de permutation issue du pivotage.
2.  **Résolution par Substitution :** Le système `P*A*x = P*b` devient `L*U*x = P*b`.
    *   On résout d'abord `L*y = P*b` pour `y` par substitution directe (forward substitution).
    *   Puis on résout `U*x = y` pour `x` par substitution arrière (backward substitution).

### Gestion de la Contrainte de Non-Négativité (`w_i >= 0`)

La version actuelle de `ResoudreSystemeOptimisation` applique une méthode simpliste :
1.  Elle résout le système KKT sans tenir compte explicitement de `w_i >= 0`.
2.  Si des poids `w_i` dans la solution sont négatifs, ils sont tronqués à `0.0`.
3.  Les poids (potentiellement tronqués) sont ensuite normalisés pour que leur somme soit égale à 1.

**Note :** Une approche plus rigoureuse pour les contraintes d'inégalité impliquerait une méthode d'ensemble actif (active-set method) ou des algorithmes de points intérieurs, qui modifient itérativement le système KKT pour satisfaire les contraintes d'inégalité. La boucle `do-while` dans `ResoudreSystemeOptimisation` était initialement prévue pour une telle méthode, mais avec `MAX_ITERATIONS_OPTIM` fixé à 1, elle ne réalise qu'une seule passe de résolution et de troncature.

## Comment Utiliser

1.  **Compiler le Code :**
    ```bash
    gcc -o optimiseur_portefeuille votre_fichier_source.c -lm
    ```
    (Remplacez `votre_fichier_source.c` par le nom de votre fichier C. `-lm` est nécessaire pour la bibliothèque mathématique).

2.  **Exécuter :**
    ```bash
    ./optimiseur_portefeuille
    ```

3.  **Personnalisation :**
    *   **Données d'Entrée :** Modifiez la section `main` pour charger vos propres données de prix d'actifs (`tempPrix`), le nombre de jours (`nombreJours`), et le nombre d'actifs (`nombreActifs`).
    *   **Rendement Cible :** Ajustez la variable `mu_target_par_periode` dans `main` pour définir le rendement attendu par période que vous visez pour le portefeuille.
    *   **Formulation KKT :** Si votre problème d'optimisation a une formulation légèrement différente pour la fonction objectif (par ex., `(1/2)w'Cw` vs `w'Cw`), ajustez la construction de la matrice KKT (notamment le facteur `2` devant la matrice de covariance).

## Applications dans d'Autres Domaines

Le noyau de ce code – la capacité à résoudre un système d'équations linéaires (`ResoudreSystemeLineaire`) et la structure pour formuler un problème d'optimisation via KKT (`ConstruireMatriceKKT`, `ResoudreSystemeOptimisation`) – peut être adapté à divers problèmes d'optimisation quadratique sous contraintes linéaires d'égalité.

Voici quelques exemples :

1.  **Ajustement de Courbes (Least Squares avec Contraintes) :**
    *   **Problème :** Trouver les paramètres d'un modèle qui minimisent la somme des carrés des erreurs (par exemple, `(Ax-b)^T(Ax-b)`), sous certaines contraintes linéaires sur les paramètres `x`.
    *   **Adaptation :** La fonction objectif est quadratique. Les contraintes linéaires peuvent être intégrées dans un système KKT similaire.

2.  **Allocation de Ressources :**
    *   **Problème :** Allouer des ressources limitées à différentes tâches pour maximiser une certaine utilité (qui peut être quadratique) ou minimiser un coût (quadratique), sous des contraintes de budget total ou de capacité.
    *   **Adaptation :** Si l'objectif est quadratique et les contraintes sont linéaires.

3.  **Ingénierie et Contrôle Optimal :**
    *   **Problème :** Déterminer les entrées de contrôle d'un système dynamique pour minimiser une fonction de coût quadratique (par exemple, l'énergie consommée ou l'erreur de suivi) sur un horizon temporel, sujet aux équations du système (contraintes linéaires).
    *   **Adaptation :** Les problèmes de type LQR (Linear Quadratic Regulator) mènent à des optimisations qui peuvent, dans certains cas, être résolues via des systèmes KKT.

4.  **Apprentissage Automatique (Machine Learning) :**
    *   **Problème :** Certains algorithmes, comme les Machines à Vecteurs de Support (SVM) avec un noyau linéaire, impliquent la résolution d'un problème d'optimisation quadratique.
    *   **Adaptation :** Bien que les SVM impliquent souvent des contraintes d'inégalité gérées par des solveurs plus spécialisés, les principes de base de la formulation KKT s'appliquent.

5.  **Reconstruction d'Images :**
    *   **Problème :** Reconstruire une image à partir de données incomplètes ou bruitées en minimisant une fonction d'erreur quadratique, potentiellement avec des contraintes de régularisation.

**Pour Adapter le Code :**

*   **Redéfinir la Fonction Objectif et les Contraintes :** Vous devrez modifier la fonction `ConstruireMatriceKKT` pour refléter la nouvelle fonction objectif quadratique et les nouvelles contraintes linéaires de votre problème.
*   **Interpréter la Solution :** Le vecteur `X_solution_iteration` contiendra les variables de décision de votre nouveau problème ainsi que les multiplicateurs de Lagrange.
*   **Gérer les Contraintes d'Inégalité :** Si votre problème comporte des contraintes d'inégalité significatives (autres que la simple non-négativité), le mécanisme de troncature actuel sera insuffisant. Il faudrait implémenter une méthode d'optimisation plus avancée (par exemple, ensemble actif itératif, points intérieurs).

Ce code fournit une base pour comprendre et implémenter des solveurs pour une classe spécifique de problèmes d'optimisation.
