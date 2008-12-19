within Modelica_Fluid.Test.TestComponents.Fittings;
model TestGenericPressureLoss
  import Modelica_Fluid;
  extends Modelica.Icons.Example;

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
    nPorts=1) 
    annotation (Placement(transformation(extent={{-40,40},{-20,60}}, rotation=0)));
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,1e5; 10,3e5]) 
    annotation (Placement(transformation(extent={{-80,40},{-60,60}}, rotation=0)));

  Modelica_Fluid.Sources.Boundary_pT ambient_p1(
    redeclare package Medium = Medium,
    T=Modelica.SIunits.Conversions.from_degC(80),
    p=200000) 
    annotation (Placement(transformation(extent={{60,40},{40,60}}, rotation=0)));
  Modelica_Fluid.Fittings.GenericPressureLoss pressureLoss(
    redeclare package Medium = Medium,
    m_flow_nominal=1,
    dp_nominal=100000) 
             annotation (Placement(transformation(extent={{0,40},{20,60}},
          rotation=0)));

  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent=
            {{-82,-90},{-62,-70}}, rotation=0)));
equation
  connect(p_table.y, ambient_a.p_in) 
                                    annotation (Line(points={{-59,50},{-52,50},
          {-52,58},{-42,58}}, color={0,0,127}));
  connect(ambient_a.ports[1], pressureLoss.port_a) 
                                           annotation (Line(points={{-20,50},{
          -10,50},{0,50}},
                color={0,127,255}));
  connect(pressureLoss.port_b, ambient_p1.ports[1]) 
                                              annotation (Line(points={{20,50},
          {40,50}}, color={0,127,255}));
end TestGenericPressureLoss;
