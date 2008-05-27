within Modelica_Fluid.Test.TestComponents.PressureLosses;
model TestSuddenExpansion
  extends Modelica.Icons.Example;
  replaceable package Medium = 
      Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater 
    constrainedby Modelica.Media.Interfaces.PartialMedium
    "Medium in all components"                        annotation (
    choicesAllMatching =                                                                            true);

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics),
    experiment(StopTime=10, NumberOfIntervals=10000),
    experimentSetupOutput,
    Documentation(info="<html>
</html>"));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambient_a(
                                                     redeclare package Medium
      = Medium,
    p=ambient.default_p_ambient,
    T=ambient.default_T_ambient,
    usePressureInput=true) 
    annotation (Placement(transformation(extent={{-40,40},{-20,60}}, rotation=0)));
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.9999e5; 10,1.0001e5]) 
    annotation (Placement(transformation(extent={{-80,40},{-60,60}}, rotation=0)));

  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (Placement(transformation(extent={{60,40},{40,60}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (Placement(transformation(extent={{60,10},{40,30}}, rotation=0)));
  Modelica_Fluid.PressureLosses.SuddenExpansion expansion1(
    redeclare package Medium = Medium,
    D_a=0.1,
    D_b=0.2,
    use_Re=false) 
             annotation (Placement(transformation(extent={{0,40},{20,60}},
          rotation=0)));
  Modelica_Fluid.PressureLosses.SuddenExpansion expansion2(
    redeclare package Medium = Medium,
    D_a=0.1,
    D_b=0.2,
    from_dp=false,
    use_Re=false) annotation (Placement(transformation(extent={{0,10},{20,30}}, 
          rotation=0)));

  inner Modelica_Fluid.Ambient ambient 
                                   annotation (Placement(transformation(extent=
            {{66,-42},{86,-22}}, rotation=0)));
equation
  connect(p_table.y, ambient_a.p_in) 
                                    annotation (Line(points={{-59,50},{-52,50},
          {-52,56},{-42,56}}, color={0,0,127}));
  connect(ambient_a.port, expansion1.port_a) 
                                           annotation (Line(points={{-20,50},{0,
          50}}, color={0,127,255}));
  connect(expansion1.port_b, ambient_p1.port) annotation (Line(points={{20,50},
          {40,50}}, color={0,127,255}));
  connect(expansion2.port_b, ambient_p2.port) annotation (Line(points={{20,20},
          {40,20}}, color={0,127,255}));
  connect(expansion2.port_a, ambient_a.port) 
                                           annotation (Line(points={{0,20},{-10,
          20},{-10,50},{-20,50}}, color={0,127,255}));
end TestSuddenExpansion;
