within Modelica_Fluid.Test.TestComponents.PressureLosses;
model TestSharpEdgedOrifice
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  // Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater
  // Modelica.Media.Water.WaterIF97_ph
  replaceable package Medium = Modelica.Media.IdealGases.SingleGases.O2 
    constrainedby Modelica.Media.Interfaces.PartialMedium
    "Medium in all components"                        annotation (
    choicesAllMatching =                                                                            true);

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics),
    experiment(StopTime=10, NumberOfIntervals=10000),
    experimentSetupOutput,
    Documentation(info="<html>
</html>"));
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambient_p(
                                                     redeclare package Medium
      = Medium,
    p=system.p_ambient,
    T=system.T_ambient,
    usePressureInput=true) 
    annotation (Placement(transformation(extent={{-40,40},{-20,60}}, rotation=0)));
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.1e5; 10,2e5]) 
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
  Modelica_Fluid.PressureLosses.SharpEdgedOrifice orifice1(
    redeclare package Medium = Medium,
    D_pipe=0.1,
    alpha=50,
    D_min=0.02,
    L=0.005,
    show_Re=true) 
             annotation (Placement(transformation(extent={{0,40},{20,60}},
          rotation=0)));
  Modelica_Fluid.PressureLosses.SharpEdgedOrifice orifice2(
    redeclare package Medium = Medium,
    D_pipe=0.1,
    alpha=50,
    from_dp=false,
    D_min=0.02,
    L=0.005) annotation (Placement(transformation(extent={{0,10},{20,30}},
          rotation=0)));

  inner Modelica_Fluid.System system 
    annotation (Placement(transformation(extent={{60,-68},{80,-48}}, rotation=0)));
equation
  connect(p_table.y, ambient_p.p_in) 
                                    annotation (Line(points={{-59,50},{-52,50},
          {-52,56},{-42,56}}, color={0,0,127}));
  connect(ambient_p.port, orifice1.port_a) 
                                         annotation (Line(points={{-20,50},{0,
          50}}, color={0,127,255}));
  connect(orifice1.port_b, ambient_p1.port) annotation (Line(points={{20,50},{
          40,50}}, color={0,127,255}));
  connect(ambient_p.port, orifice2.port_a) 
                                         annotation (Line(points={{-20,50},{-12,
          50},{-12,20},{0,20}}, color={0,127,255}));
  connect(orifice2.port_b, ambient_p2.port) annotation (Line(points={{20,20},{
          40,20}}, color={0,127,255}));
end TestSharpEdgedOrifice;
