within Modelica_Fluid.Test.TestComponents.PressureLosses;
model TestSharpEdgedOrifice 
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  // Modelica.Media.Water.ConstantPropertyLiquidWater 
  // Modelica.Media.Water.WaterIF97_ph 
  replaceable package Medium = Modelica.Media.IdealGases.SingleGases.O2 
    extends Modelica.Media.Interfaces.PartialMedium "Medium in all components" 
                                                      annotation (
    choicesAllMatching =                                                                            true);
  
  annotation (Diagram,
    experiment(StopTime=10, NumberOfIntervals=10000),
    experimentSetupOutput,
    Coordsys(extent=[-100,-100; 100,100]),
    Documentation(info="<html>
</html>"));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambient_p(
                                                     redeclare package Medium 
      = Medium,
    p=ambient.default_p_ambient,
    T=ambient.default_T_ambient) 
    annotation (extent=[-40,40; -20,60]);
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.1e5; 10,2e5]) 
    annotation (extent=[-80,40; -60,60]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,40; 40,60]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[60,10; 40,30]);
  Modelica_Fluid.PressureLosses.SharpEdgedOrifice orifice1(
    redeclare package Medium = Medium,
    D_pipe=0.1,
    alpha=50,
    D_min=0.02,
    L=0.005,
    show_Re=true) 
             annotation (extent=[0,40; 20,60]);
  Modelica_Fluid.PressureLosses.SharpEdgedOrifice orifice2(
    redeclare package Medium = Medium,
    D_pipe=0.1,
    alpha=50,
    from_dp=false,
    D_min=0.02,
    L=0.005) annotation (extent=[0,10; 20,30]);
  
  inner Modelica_Fluid.Ambient ambient 
    annotation (extent=[60,-68; 80,-48]);
equation 
  connect(p_table.y, ambient_p.p_in) 
                                    annotation (points=[-59,50; -52,50; -52,56;
        -42,56], style(color=74, rgbcolor={0,0,127}));
  connect(ambient_p.port, orifice1.port_a) 
                                         annotation (points=[-20,50; 0,50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
  connect(orifice1.port_b, ambient_p1.port) annotation (points=[20,50; 40,
        50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
  connect(ambient_p.port, orifice2.port_a) 
                                         annotation (points=[-20,50; -12,
        50; -12,20; 0,20],
                        style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
  connect(orifice2.port_b, ambient_p2.port) annotation (points=[20,20; 40,
        20],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=3,
      rgbfillColor={0,0,255},
      fillPattern=8));
end TestSharpEdgedOrifice;
