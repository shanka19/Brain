# Debug Levels
#  bits: xxxx
#  bit 0 = TCanvas visualization
#  bit 1 = Verbalization
#  bit 2 = Fill histograms
#  bit 3 = Draw the histograms

./BrainInWorld -debug 0 -skipGenerations 500 -endGeneration 500 -timeStep 200 -worldSize 500 -nBots 10 -nFoods 5 -mu_modConnection 0.05 -mu_visualAngle 0.008 -reproductionType 2
