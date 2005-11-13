model TestSharpEdgedOrifice 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  replaceable package Medium = 
      Modelica.Media.Water.ConstantPropertyLiquidWater 
    extends Modelica.Media.Interfaces.PartialMedium "Medium in all components" 
                                                      annotation (
    choicesAllMatching =                                                                            true);
  
  annotation (Diagram,
    experiment(StopTime=10, NumberOfIntervals=10000),
    experimentSetupOutput,
    Coordsys(extent=[-100,-100; 100,100]),
    Documentation(info="<html>
</html>"));
  Sources.PrescribedAmbient_pTX ambient(redeclare package Medium = Medium) 
    annotation (extent=[-40,40; -20,60]);
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.1e5; 10,2e5]) 
    annotation (extent=[-80,40; -60,60]);
  
  Sources.FixedAmbient_pTX ambient_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[62,40; 42,60]);
  Sources.FixedAmbient_pTX ambient_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,10; 40,30]);
  Modelica_Fluid.WorkInProgress.Components.SharpEdgedOrifice orifice1(
    redeclare package Medium = Medium,
    D_pipe=0.1,
    alpha=50,
    D_min=0.02,
    L=0.005) annotation (extent=[0,40; 20,60]);
  Modelica_Fluid.WorkInProgress.Components.SharpEdgedOrifice orifice2(
    redeclare package Medium = Medium,
    D_pipe=0.1,
    alpha=50,
    from_dp=false,
    D_min=0.02,
    L=0.005) annotation (extent=[0,10; 20,30]);
equation 
  connect(p_table.y, ambient.p_in)  annotation (points=[-59,50; -52,50; -52,56;
        -42,56], style(color=74, rgbcolor={0,0,127}));
  connect(ambient.port, orifice1.port_a) annotation (points=[-19,50; -1,50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
  connect(orifice1.port_b, ambient_p1.port) annotation (points=[21,50; 41,50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
  connect(ambient.port, orifice2.port_a) annotation (points=[-19,50; -12,50;
        -12,20; -1,20], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
  connect(orifice2.port_b, ambient_p2.port) annotation (points=[21,20; 39,20],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
end TestSharpEdgedOrifice;
