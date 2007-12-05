within Modelica_Fluid.Test.TestComponents.PressureLosses;
model TestSuddenExpansion 
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
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambient_a(
                                                     redeclare package Medium 
      = Medium,
    p=ambient.default_p_ambient,
    T=ambient.default_T_ambient) 
    annotation (extent=[-40,40; -20,60]);
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.9999e5; 10,1.0001e5]) 
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
  Modelica_Fluid.PressureLosses.SuddenExpansion expansion1(
    redeclare package Medium = Medium,
    D_a=0.1,
    D_b=0.2,
    use_Re=false) 
             annotation (extent=[0,40; 20,60]);
  Modelica_Fluid.PressureLosses.SuddenExpansion expansion2(
    redeclare package Medium = Medium,
    D_a=0.1,
    D_b=0.2,
    from_dp=false,
    use_Re=false) annotation (extent=[0,10; 20,30]);
  
  inner Modelica_Fluid.Ambient ambient 
                                   annotation (extent=[66,-42; 86,-22]);
equation 
  connect(p_table.y, ambient_a.p_in) 
                                    annotation (points=[-59,50; -52,50; -52,56;
        -42,56], style(color=74, rgbcolor={0,0,127}));
  connect(ambient_a.port, expansion1.port_a) 
                                           annotation (points=[-20,50; 0,50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
  connect(expansion1.port_b, ambient_p1.port) annotation (points=[20,50; 40,50],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
  connect(expansion2.port_b, ambient_p2.port) annotation (points=[20,20; 40,20],
      style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
  connect(expansion2.port_a, ambient_a.port) 
                                           annotation (points=[0,20; -10,20;
        -10,50; -20,50], style(
      color=69,
      rgbcolor={0,127,255},
      fillColor=7,
      rgbfillColor={255,255,255},
      fillPattern=1));
end TestSuddenExpansion;
