within Modelica_Fluid.Test.TestComponents.Fittings;
model TestSimpleGenericOrifice
  extends Modelica.Icons.Example;
  parameter Real diameter_a = 0.1;
  parameter Real diameter_b = 0.2;
  parameter Real A_rel = (diameter_a/diameter_b)^2;
  parameter Real zeta =  (1 - A_rel)^2;

  replaceable package Medium = 
      Modelica.Media.Water.ConstantPropertyLiquidWater 
    constrainedby Modelica.Media.Interfaces.PartialMedium
    "Medium in all components"                        annotation (
    choicesAllMatching =                                                                            true);

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics),
    experiment(StopTime=10, NumberOfIntervals=10000),
    experimentSetupOutput,
    Documentation(info="<html>
</html>"));
  Modelica_Fluid.Sources.Boundary_pT ambient_a(     redeclare package Medium = 
        Medium,
    p=system.p_ambient,
    T=system.T_ambient,
    use_p_in=true,
    nPorts=4) 
    annotation (Placement(transformation(extent={{-40,40},{-20,60}}, rotation=0)));
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.9999e5; 10,1.0001e5]) 
    annotation (Placement(transformation(extent={{-80,40},{-60,60}}, rotation=0)));

  Modelica_Fluid.Sources.Boundary_pT ambient_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (Placement(transformation(extent={{60,40},{40,60}}, rotation=0)));
  Modelica_Fluid.Sources.Boundary_pT ambient_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (Placement(transformation(extent={{60,10},{40,30}}, rotation=0)));
  Modelica_Fluid.Fittings.SimpleGenericOrifice expansion1(
    redeclare package Medium = Medium,
    zeta=zeta,
    diameter=diameter_a) 
             annotation (Placement(transformation(extent={{0,40},{20,60}},
          rotation=0)));
  Modelica_Fluid.Fittings.SimpleGenericOrifice expansion2(
    redeclare package Medium = Medium,
    from_dp=false,
    zeta=zeta,
    diameter=diameter_a) annotation (Placement(transformation(extent={{0,10},{20,30}},
          rotation=0)));

  Modelica_Fluid.Fittings.SuddenExpansion expansion3(
    redeclare package Medium = Medium,
    use_Re=false,
    diameter_a=diameter_a,
    diameter_b=diameter_b) annotation (Placement(transformation(extent={{0,-20},{20,0}},
          rotation=0)));
  Modelica_Fluid.Sources.Boundary_pT ambient_p3(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (Placement(transformation(extent={{60,-20},{40,0}}, rotation=0)));
  Modelica_Fluid.Fittings.SuddenExpansion expansion4(
    redeclare package Medium = Medium,
    use_Re=false,
    diameter_a=diameter_a,
    diameter_b=diameter_b,
    from_dp=false) 
             annotation (Placement(transformation(extent={{0,-60},{20,-40}},
          rotation=0)));
  Modelica_Fluid.Sources.Boundary_pT ambient_p4(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (Placement(transformation(extent={{60,-60},{40,-40}}, rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent=
            {{-82,-90},{-62,-70}}, rotation=0)));
equation
  connect(p_table.y, ambient_a.p_in) 
                                    annotation (Line(points={{-59,50},{-52,50},
          {-52,58},{-42,58}}, color={0,0,127}));
  connect(ambient_a.ports[1], expansion1.port_a) 
                                           annotation (Line(points={{-20,53},{
          -10,53},{-10,50},{0,50}},
                color={0,127,255}));
  connect(expansion1.port_b, ambient_p1.ports[1]) 
                                              annotation (Line(points={{20,50},
          {40,50}}, color={0,127,255}));
  connect(expansion2.port_b, ambient_p2.ports[1]) 
                                              annotation (Line(points={{20,20},
          {40,20}}, color={0,127,255}));
  connect(expansion2.port_a, ambient_a.ports[2]) 
                                           annotation (Line(points={{0,20},{-10,
          20},{-10,51},{-20,51}}, color={0,127,255}));
  connect(ambient_a.ports[3], expansion3.port_a) 
                                           annotation (Line(points={{-20,49},{
          -10,49},{-10,-10},{0,-10}}, color={0,127,255}));
  connect(expansion3.port_b, ambient_p3.ports[1]) 
                                              annotation (Line(points={{20,-10},
          {40,-10}}, color={0,127,255}));
  connect(ambient_a.ports[4], expansion4.port_a) 
                                           annotation (Line(points={{-20,47},{
          -10,47},{-10,-50},{0,-50}}, color={0,127,255}));
  connect(expansion4.port_b, ambient_p4.ports[1]) 
                                              annotation (Line(points={{20,-50},
          {40,-50}}, color={0,127,255}));
end TestSimpleGenericOrifice;
