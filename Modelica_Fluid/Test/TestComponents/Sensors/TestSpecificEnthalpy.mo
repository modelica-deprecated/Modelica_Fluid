within Modelica_Fluid.Test.TestComponents.Sensors;
model TestSpecificEnthalpy 
  import Modelica_Fluid;
  annotation (
    Diagram,
    experiment(Tolerance=1e-006),
    experimentSetupOutput);
  inner Modelica_Fluid.Ambient ambient annotation (extent=[-100,-100; -80,-80]);
  Modelica_Fluid.Sources.PrescribedBoundary_phX boundary_prescribed_1(
    useEnthalpyInput=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient) annotation (extent=[-40,10; -20,30]);
  Modelica_Fluid.Sensors.SpecificEnthalpyOnePort specificEnthalpy(redeclare 
      package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-10,20; 10,40]);
  Modelica_Fluid.Sources.PrescribedBoundary_phX boundary_prescribed_2(
    useEnthalpyInput=true,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=ambient.default_p_ambient) annotation (extent=[-40,-30; -20,-10]);
  Modelica_Fluid.Sensors.SpecificEnthalpyTwoPort specificEnthalpy1(redeclare 
      package Medium = Modelica.Media.Water.StandardWater) 
    annotation (extent=[-10,-30; 10,-10]);
  Modelica.Blocks.Sources.Sine sine1(amplitude=1600e3, offset=1800e3) 
                                    annotation (extent=[-80,-10; -60,10]);
equation 
  connect(boundary_prescribed_1.port, specificEnthalpy.port) annotation (points=
       [-20,20; 0,20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(sine1.y, boundary_prescribed_1.h_in) annotation (points=[-59,0; -52,0;
        -52,20; -42,20], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(sine1.y, boundary_prescribed_2.h_in) annotation (points=[-59,0; -52,0;
        -52,-20; -42,-20], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(boundary_prescribed_2.port, specificEnthalpy1.port_a) annotation (
      points=[-20,-20; -10,-20], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
end TestSpecificEnthalpy;
