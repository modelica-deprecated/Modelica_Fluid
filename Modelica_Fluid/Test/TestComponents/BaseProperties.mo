within Modelica_Fluid.Test.TestComponents;
package BaseProperties "Test BaseProperties models"
  model TestBaseProperties01 "Simple test of a single BaseProperties model"

    replaceable package Medium = 
          ModelicaNew.Media.FluidTypes.IdealGases.NasaAirN2O2 constrainedby
      ModelicaNew.Media.Interfaces.GenericMedium annotation (choicesAllMatching = true);

    Media.BaseProperties medium(redeclare package Medium = 
          Medium) 
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  equation
    medium.p = Medium.p_default;
    medium.T = Medium.T_default*(1 + 0.01*Modelica.Math.exp(time/10)*Modelica.Math.sin(2*3.14*time));
    medium.X = Medium.X_default;

    annotation (experiment(StopTime=10));
  end TestBaseProperties01;



end BaseProperties;
