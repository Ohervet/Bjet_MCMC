{

gROOT->Reset();

//gSystem->Exec("./bj PKS_0625.par");
//gSystem->Exec("./bj BLLAC.par");
//gSystem->Exec("./bj W_Comae.par");
//gSystem->Exec("./bj VER_J0521.par");
//gSystem->Exec("./bj HESS_J1943.par");
//gSystem->Exec("./bj PMN_J1603.par");
//gSystem->Exec("./bj 1ES_1215.par");
gSystem->Exec("./bj OJ_287.par");
//gSystem->Exec("./bj NGC_1275.par");
//gSystem->Exec("./bj TON_599.par");
//gSystem->Exec("./bj PKS_1222.par");
//gSystem->Exec("./gal galaxy.par");
gROOT->ProcessLine(".x plmodel.C");
}
