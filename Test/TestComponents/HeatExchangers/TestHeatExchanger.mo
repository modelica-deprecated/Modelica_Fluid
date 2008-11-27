within Modelica_Fluid.Test.TestComponents.HeatExchangers;
model TestHeatExchanger

extends Modelica.Icons.Example;

//replaceable package Medium = Modelica.Media.Water.StandardWater;
package Medium = Modelica.Media.Incompressible.Examples.Essotherm650;
  Modelica_Fluid.HeatExchangers.BasicHX HEX(
    c_wall=500,
    use_T_start=true,
    n=20,
    length=2,
    m_flow_start_1=0.2,
    m_flow_start_2=0.2,
    redeclare package WallFriction_1 = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    redeclare package WallFriction_2 = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
    k_wall=100,
    initType=Modelica_Fluid.Types.Init.SteadyStateHydraulic,
    s_wall=0.005,
    Ac_1=4.5e-4,
    Ac_2=4.5e-4,
    P_1=0.075,
    P_2=0.075,
    Ah_1=0.075*2*20,
    Ah_2=0.075*2*20,
    d_wall=900,
    redeclare package Medium_1 = 
        Medium,
    redeclare package Medium_2 = 
        Medium,
    redeclare model HeatTransfer_1 = 
        Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha (alpha0=
           1000),
    Twall_start=300,
    dT=10,
    T_start_1=304,
    T_start_2=300,
    eta_nominal_1=0.01,
    eta_nominal_2=0.01,
    use_eta_nominal=true)       annotation (Placement(transformation(extent={{
            -26,-14},{34,46}}, rotation=0)));

  Modelica_Fluid.Sources.FixedBoundary_pTX ambient2(
    p=1e5,
    T=280,
    redeclare package Medium = Medium)                              annotation (Placement(
        transformation(extent={{82,-28},{62,-8}}, rotation=0)));
  Modelica_Fluid.Sources.FixedBoundary_pTX ambient1(
    p=1e5,
    T=300,
    redeclare package Medium = Medium)                              annotation (Placement(
        transformation(extent={{82,24},{62,44}}, rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate2(
    m_flow=0.2,
    T=360,
    redeclare package Medium = Medium,
    useFlowRateInput=true,
    useTemperatureInput=false,
    useCompositionInput=false) 
                annotation (Placement(transformation(extent={{-66,24},{-46,44}},
          rotation=0)));
  Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowRate1(
    T=300,
    m_flow=0.5,
    redeclare package Medium = Medium) 
                 annotation (Placement(transformation(extent={{-66,-10},{-46,10}},
          rotation=0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}),
                      graphics),
                       experiment(StopTime=100, Tolerance=1e-005));
  Modelica.Blocks.Sources.Ramp Ramp1(
    startTime=50,
    duration=5,
    height=-1,
    offset=0.5)   annotation (Placement(transformation(extent={{-100,24},{-80,
            44}}, rotation=0)));
  inner Modelica_Fluid.System system 
                                   annotation (Placement(transformation(extent=
            {{60,70},{80,90}}, rotation=0)));
equation
  connect(massFlowRate1.port, HEX.port_a1)            annotation (Line(points={
          {-46,0},{-40,0},{-40,15.4},{-29,15.4}}, color={0,127,255}));
  connect(HEX.port_b1, ambient1.port)            annotation (Line(points={{37,
          15.4},{48.5,15.4},{48.5,34},{62,34}}, color={0,127,255}));
  connect(Ramp1.y, massFlowRate2.m_flow_in) annotation (Line(points={{-79,34},{
          -74,34},{-74,40},{-65.3,40}}, color={0,0,127}));
  connect(massFlowRate2.port, HEX.port_b2) annotation (Line(
      points={{-46,34},{-40,34},{-40,29.8},{-29,29.8}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(HEX.port_a2, ambient2.port) annotation (Line(
      points={{37,2.2},{42,2},{50,2},{50,-18},{62,-18}},
      color={0,127,255},
      smooth=Smooth.None));
end TestHeatExchanger;
