
exec("./scilab-scripts/plotRIMsolutions.sce");

//Plot benchmark
scf(0);clf();
plotRIM_K_L();

scf(1);clf();
plotRIM_Kir_K_L();

scf(2);clf();
plotRIM_Cap_K_L();

scf(3);clf();
plotRIM_Cap_Kir_K_L();

//Plot benchmark with found best solutions
scf(4);clf();
plotRIM_K_L(%T);

scf(5);clf();
plotRIM_Kir_K_L(%T);

scf(6);clf();
plotRIM_Cap_K_L(%T);

scf(7);clf();
plotRIM_Cap_Kir_K_L(%T);

//Export benchmark figures to pdf

xs2pdf(0,'RIM_K_L.pdf');
xs2pdf(1,'RIM_Kir_K_L.pdf');
xs2pdf(2,'RIM_Cap_K_L.pdf');
xs2pdf(3,'RIM_Cap_Kir_K_L.pdf');

//Export results figures to pdf

xs2pdf(4,'resultsRIM_K_L.pdf');
xs2pdf(5,'resultsRIM_Kir_K_L.pdf');
xs2pdf(6,'resultsRIM_Cap_K_L.pdf');
xs2pdf(7,'resultsRIM_Cap_Kir_K_L.pdf');
