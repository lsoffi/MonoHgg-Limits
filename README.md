Instructions on how to run the full limit computation and plotting

----------------------------------------------------------
#2HDM
 ./combineall_MonoHgg.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --suffix ""
ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9

 ./combineall_MonoHgg.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --suffix "_HighMET"

 ./combineall_MonoHgg.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --suffix "_LowMET"

Make all limit plots:
source limit_plots_MonoHgg.sh 

----------------------------------------------------------
#Baryonic
./combineall_MonoHgg_ZpBaryonic.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --suffix ""

./combineall_MonoHgg_ZpBaryonic.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --suffix "_HighMET"

./combineall_MonoHgg_ZpBaryonic.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --suffix "_LowMET"

Make all limit plots:
source limit_plots_MonoHgg_ZpBaryonic.sh 

----------------------------------------------------------

#Scalar

./combineall_MonoHgg_Scalar.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --suffix ""

./combineall_MonoHgg_Scalar.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --_

./combineall_MonoHgg_Scalar.sh ntuples4fit_pho_met0_met130_cic_default_shapes_lumi_35.9 --hadd -M Asymptotic --run both --_HighMET"

Make all limit plots:
source limit_plots_MonoHgg_Scalar.sh 

----------------------------------------

