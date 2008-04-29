within Modelica_Fluid.Test.TestComponents.Utilities;
model TestRegStep "Test regStep and regUnitStep functions" 
  extends Modelica.Icons.Example;
  import Modelica_Fluid.Utilities.*;
  parameter Real x_small=0.5;
  Real x=time - 0.6;
  Real yRegUnitStep = Modelica_Fluid.Utilities.regUnitStep(x,x_small);
  Real yRegStep = Modelica_Fluid.Utilities.regStep(x,1,0.5,x_small);
  Real yRegStep_der;
  annotation (experiment(StopTime=1.3), experimentSetupOutput);
equation 
  yRegStep_der = der(yRegStep);
end TestRegStep;
