// Création de la frontière du triangle
border b1(t=0,1){x=t; y=0;};    // Segment de (0,0) à (1,0)
border b2(t=0,1){x=1-t; y=t;};  // Segment de (1,0) à (0,1)
border b3(t=0,1){x=0; y=1-t;};  // Segment de (0,1) à (0,0)

// Définition des solutions exactes
func uexact = x;
func vexact = -y;
func pexact = x + y + 1;

// Initialisation des forces extérieures
real f1 = 0.0, f2 = 0.0;

// Nombre d'itérations de raffinement
int numRaffinements = 5;
real[int] erreursU(numRaffinements);
real[int] erreursP(numRaffinements);
int[int] nbTriangles(numRaffinements);

// Boucle pour raffiner le maillage et calculer l'erreur
for (int i = 0; i < numRaffinements; ++i) {
    // Création du maillage avec raffinement
    mesh Th = buildmesh(b1(10*(i+1)) + b2(10*(i+1)) + b3(10*(i+1)));
    nbTriangles[i] = Th.nt;  // Enregistrer le nombre de triangles

    // Espaces finis
    fespace Uh(Th, P1b); // Espace pour la vitesse (P1 bulle)
    fespace Ph(Th, P1);  // Espace pour la pression (P1)

    Uh u, v, uu, vv;  // Composantes de la vitesse
    Ph p, pp;         // Pression

    // Problème de Stokes
    solve stokes([u, v, p], [uu, vv, pp], solver = UMFPACK)
        = int2d(Th)(
              dx(u)*dx(uu) + dy(u)*dy(uu)
            + dx(v)*dx(vv) + dy(v)*dy(vv)
            - dx(uu)*p - dy(vv)*p
            + pp*(dx(u) + dy(v))
            - 1e-6*p*pp
          )
        + int2d(Th)(f1*uu + f2*vv)
        + on(1, u=uexact, v=vexact)
        + on(2, u=uexact, v=vexact)
        + on(3, u=uexact, v=vexact);

    // Calcul de l'erreur en norme H1 pour la vitesse
    erreursU[i] = sqrt(int2d(Th)(
        (u - uexact)^2 + (v - vexact)^2
        + (dx(u) - 1)^2 + (dy(u) - 0)^2
        + (dx(v) - 0)^2 + (dy(v) + 1)^2
    ));

    // Calcul de l'erreur en L2 pour la pression
    erreursP[i] = sqrt(int2d(Th)((p - pexact)^2));
}

// Affichage des erreurs
for (int i = 0; i < numRaffinements; ++i) {
    cout << "Raffinement " << i+1 << ", Nombre de triangles: " << nbTriangles[i]
         << ", Erreur u: " << erreursU[i]
         << ", Erreur p: " << erreursP[i] << endl;
}

// Graphe des erreurs
plot(nbTriangles, erreursU, cmm="Erreur sur u en fonction du nombre de triangles", wait=false);
plot(nbTriangles, erreursP, cmm="Erreur sur p en fonction du nombre de triangles");


