#include "GeneratorPV.hh"
//#include "Generator.hh"
#include "PhysicalConstants.hh" // for HBARC and M_NUCLEON
//#include "Commutator.hh"
//#include "Operator.hh"
//GeneratorPV::GeneratorPV():
//             Generator() 
//{}



using PhysConst::M_NUCLEON;
using PhysConst::HBARC;



void GeneratorPV::Update(Operator& H_s, Operator& HPV_s, Operator& Eta_s, Operator& EtaPV_s)
{   
   Eta_s.Erase();
   EtaPV_s.Erase();
   AddToEtaPV(H_s,HPV_s,Eta_s,EtaPV_s);

}

void GeneratorPV::AddToEtaPV(Operator& H_s, Operator& HPV_s, Operator& Eta_s, Operator& EtaPV_s)
{
   double start_time = omp_get_wtime();
   H = &H_s;
   Eta = &Eta_s;
   V = &HPV_s;
   Etapv = &EtaPV_s;
   modelspace = H->GetModelSpace();

        if (generator_type == "wegner")                       ConstructGeneratorPV_SingleRef(wegner_func); // never tested, probably doesn't work.
   else if (generator_type == "white")                        ConstructGeneratorPV_SingleRef(white_func);
   else if (generator_type == "atan")                         ConstructGeneratorPV_SingleRef(atan_func) ;
   else if (generator_type == "imaginary-time")               ConstructGeneratorPV_SingleRef(imaginarytime_func);
   else if (generator_type == "qtransfer-atan")               ConstructGeneratorPV_SingleRef(qtransferatan1_func);
//   else if (generator_type == "shell-model-wegner")           ConstructGeneratorPV_ShellModel(wegner_func); // never tested, probably doesn't work.
//   else if (generator_type == "shell-model")                  ConstructGeneratorPV_ShellModel(white_func);
//   else if (generator_type == "shell-model-atan")             ConstructGeneratorPV_ShellModel(atan_func);
//   else if (generator_type == "shell-model-atan-npnh")        ConstructGeneratorPV_ShellModel_NpNh(atan_func);
//   else if (generator_type == "shell-model-imaginary-time")   ConstructGeneratorPV_ShellModel(imaginarytime_func);
//   else if (generator_type == "hartree-fock")                 ConstructGeneratorPV_HartreeFock();
//   else if (generator_type == "1PA")                          ConstructGeneratorPV_1PA(atan_func);
   else if (generator_type.find("qtransfer-atan") != std::string::npos )
   {
     int n;
     std::istringstream( generator_type.substr( generator_type.find("_")+1) ) >> n;
     std::function<double(double,double)> qtransferatanN_func = [n](double Hod, double denom){return pow(std::abs(denom)*M_NUCLEON/HBARC/HBARC, 0.5*n) * atan_func(Hod, denom);};
//     ConstructGenerator_QTransferAtan(n);
     ConstructGeneratorPV_SingleRef( qtransferatanN_func );
   }
//   else if (generator_type == "rspace")                       ConstructGenerator_Rspace();
   else
   {
      std::cout << "Error. Unkown generator_type: " << generator_type << std::endl;
   }
   Eta->profiler.timer["UpdateEta"] += omp_get_wtime() - start_time;
   Etapv->profiler.timer["UpdateEta"] += omp_get_wtime();
}






void GeneratorPV::ConstructGeneratorPV_SingleRef(std::function<double (double,double)>& etafunc )
{
      // One body piece -- eliminate ph bits
   for ( auto& a : modelspace->core)
   {
      for ( auto& i : VectorUnion(modelspace->valence,modelspace->qspace) )
      {
         double denominator = Get1bDenominator(i,a);
//         Eta->OneBody(i,a) = 0.5*atan(2*H->OneBody(i,a)/denominator);
//         std::cout<<"denominator="<<denominator<<__LINE__<<std::endl;
         Eta->OneBody(i,a) = etafunc( H->OneBody(i,a), denominator);
         Eta->OneBody(a,i) = - Eta->OneBody(i,a);
         Etapv->OneBody(i,a) = etafunc( V->OneBody(i,a), denominator);
         Etapv->OneBody(a,i) = - Etapv->OneBody(i,a); //old version
//         std::cout<<"numeratorH="<<H->OneBody(i,a)<<__LINE__<<std::endl;
//         std::cout<<"i orbit="<< i << " " <<__LINE__<<std::endl;
//         std::cout<<"a orbit="<< a << " " <<__LINE__<<std::endl;
//         std::cout<<"numeratorV="<<V->OneBody(i,a)<<__LINE__<<std::endl;
//         std::cout<<"i orbit="<< i << " " <<__LINE__<<std::endl;
//         std::cout<<"a orbit="<< a << " " <<__LINE__<<std::endl;
     }
   }
   for ( auto& iter : Eta->TwoBody.MatEl )
  {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
//      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
//      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& ETA2 =  iter.second;
//      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
    arma::mat& H2 = H->TwoBody.GetMatrix(ch_bra,ch_ket);
//      for ( auto& iket : tbc.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)

    for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
//         for ( auto& ibra : VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
         {
//            double denominator = Get2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator(ch_bra,ch_ket,ibra,iket);
//            double denominator = Get2bDenominator_Jdep(ch,ibra,iket);
//            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETA2(ibra,iket) = etafunc( H2(ibra,iket), denominator);
            ETA2(iket,ibra) = - ETA2(ibra,iket) ; // Eta needs to be antisymmetric
         }
      }
    }
   for ( auto& iter : Etapv->TwoBody.MatEl )
   {
      size_t ch_bra = iter.first[0];
      size_t ch_ket = iter.first[1];
//      TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
      TwoBodyChannel& tbc_bra = modelspace->GetTwoBodyChannel(ch_bra);
      TwoBodyChannel& tbc_ket = modelspace->GetTwoBodyChannel(ch_ket);
//      arma::mat& ETA2 =  Eta->TwoBody.GetMatrix(ch);
      arma::mat& ETAPV2 =  iter.second;
//      arma::mat& H2 = H->TwoBody.GetMatrix(ch);
    arma::mat& V2 = V->TwoBody.GetMatrix(ch_bra,ch_ket);
//      for ( auto& iket : tbc.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)

    for ( auto& iket : tbc_ket.GetKetIndex_cc() ) // cc means core-core ('holes' refer to the reference state)
      {
//         for ( auto& ibra : VectorUnion(tbc.GetKetIndex_qq(), tbc.GetKetIndex_vv(), tbc.GetKetIndex_qv() ) )
         for ( auto& ibra : VectorUnion(tbc_bra.GetKetIndex_qq(), tbc_bra.GetKetIndex_vv(), tbc_bra.GetKetIndex_qv() ) )
             { 
//            double denominator = Get2bDenominator(ch,ibra,iket);
            double denominator = Get2bDenominator(ch_bra,ch_ket,ibra,iket);
//            double denominator = Get2bDenominator_Jdep(ch,ibra,iket);
//            ETA2(ibra,iket) = 0.5*atan(2*H2(ibra,iket) / denominator);
            ETAPV2(ibra,iket) = etafunc( V2(ibra,iket), denominator);
 //   std::cout<<__LINE__<<std::endl;
//           ETAPV2(iket,ibra) = - ETAPV2(ibra,iket) ; // Eta needs to be antisymmetric error comes from here
//           ETAPV2(iket,ibra) = etafunc( V2(iket,ibra), denominator);
//            std::cout<<"twobodyV"<< V2 << std::endl;

        }
      }
    }


}










