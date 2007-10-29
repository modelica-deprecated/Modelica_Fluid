within Modelica_Fluid.Test.TestComponents.Volumes;
model TestSweptVolume 
  "Enclosed medium with fixed quantity in an adiabatic volume with varying size" 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  Modelica_Fluid.Volumes.SweptVolume sweptVolume(
    pistonCrossArea=0.05*0.05*Modelica.Constants.pi/4, 
    clearance=0.05*0.05*Modelica.Constants.pi/4*0.03, 
    redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
    annotation (extent=[0,-10; 20,10]);
  Modelica.Mechanics.Translational.Position position 
    annotation (extent=[-20,-40; 0,-20]);
  Modelica.Blocks.Sources.Sine sine(
    phase=Modelica.Constants.pi/2,
    amplitude=0.1,
    offset=0.1) annotation (extent=[-80,-40; -60,-20]);
  inner Modelica_Fluid.Ambient ambient annotation (extent=[80,-40; 100,-20]);
  
  annotation (Diagram(Text(
        extent=[-100,80; 100,60], 
        style(color=0, rgbcolor={0,0,0}), 
        string=
            "Enclosed medium with fixed quantity in an adiabatic volume with varying size")), 
      
    experiment(StopTime=10, Tolerance=1e-007), 
    experimentSetupOutput);
equation 
  connect(position.flange_b, sweptVolume.flange) annotation (points=[0,-30; 10,
        -30; 10,-10],            style(
      color=58,
      rgbcolor={0,127,0},
      pattern=0,
      smooth=0,
      fillColor=9,
      rgbfillColor={135,135,135},
      fillPattern=7));
  connect(sine.y, position.s_ref) annotation (points=[-59,-30; -22,-30], style(
      color=74,
      rgbcolor={0,0,127},
      pattern=0,
      smooth=0,
      fillColor=9,
      rgbfillColor={135,135,135},
      fillPattern=7));
  
end TestSweptVolume;
