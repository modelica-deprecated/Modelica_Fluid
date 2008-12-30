within Modelica_Fluid.Examples;
model InverseParameterization
  "Demonstrates the parameterization of a pump and a pipe for given nominal values"
  import Modelica_Fluid;
  extends Modelica.Icons.Example;
  replaceable package Medium = Modelica.Media.Water.StandardWater;
      //Modelica.Media.Water.ConstantPropertyLiquidWater;

  Modelica_Fluid.Sources.Boundary_pT source(
    redeclare package Medium = Medium,
    nPorts=1,
    use_p_in=false,
    p=100000) 
    annotation (Placement(transformation(extent={{-76,-6},{-64,6}},  rotation=0)));
  Modelica_Fluid.Machines.ControlledPump pump(
    m_flow_nominal=1,
    control_m_flow=false,
    use_p_set=true,
    redeclare package Medium = Medium,
    p_a_nominal=100000,
    p_b_nominal=200000) 
    annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
  Modelica_Fluid.Fittings.SimpleGenericOrifice orifice(
    redeclare package Medium = Medium,
    diameter=2.54e-2,
    m_flow_nominal=1,
    use_zeta=false,
    zeta=0,
    dp_nominal=100000) 
                      annotation (Placement(transformation(extent={{20,-10},{40,
            10}}, rotation=0)));

  Modelica_Fluid.Sources.Boundary_pT sink1(redeclare package Medium = Medium, p=
       100000) 
             annotation (Placement(transformation(extent={{76,-6},{64,6}},
          rotation=0)));

  inner Modelica_Fluid.System system 
                        annotation (Placement(transformation(extent={{-90,70},{
            -70,90}},  rotation=0)));
  Modelica_Fluid.Pipes.StaticPipe pipe(
    redeclare package Medium = Medium,
    diameter=2.54e-2,
    length=0,
    redeclare model PressureLoss = 
        Modelica_Fluid.Pipes.BaseClasses.PressureLoss.NominalTurbulentFlow (
        dp_nominal=100000,
        m_flow_nominal=1,
        show_Res=true)) 
                      annotation (Placement(transformation(extent={{20,-50},{40,
            -30}},rotation=0)));
  Modelica_Fluid.Sources.Boundary_pT sink2(
    redeclare package Medium = Medium, p=200000) 
             annotation (Placement(transformation(extent={{76,-46},{64,-34}},
          rotation=0)));
  Modelica.Blocks.Sources.Ramp p_set(
    height=0.2e5,
    offset=1.9e5,
    duration=8,
    startTime=1) 
    annotation (Placement(transformation(extent={{-60,30},{-40,50}})));
equation
  connect(orifice.port_b, sink1.ports[1]) 
                                       annotation (Line(
      points={{40,0},{64,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(source.ports[1], pump.port_a) annotation (Line(
      points={{-64,0},{-40,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pump.port_b, orifice.port_a) 
                                    annotation (Line(
      points={{-20,0},{20,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe.port_b, sink2.ports[1])  annotation (Line(
      points={{40,-40},{64,-40}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(pipe.port_a, pump.port_b)  annotation (Line(
      points={{20,-40},{0,-40},{0,0},{-20,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(p_set.y, pump.p_set) annotation (Line(
      points={{-39,40},{-25,40},{-25,8.2}},
      color={0,0,127},
      smooth=Smooth.None));

  annotation (
    Commands(file(ensureSimulated=true)="Scripts/Examples/InverseParameterization/plotResults.mos"
        "plotResults"),
    Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
            100,100}}),
            graphics),
    experiment(StopTime=10, NumberOfIntervals=10000),
    Documentation(info="<html>
<p>
A pump, an orifice and a pipe are parameterized with simple nominal values. 
Note that the pipe uses the pressureLoss model NominalTurbulentFlow, which does not require the specification of geometry data. 
Instead it internally parameterizes a QuadraticTurbulentFlow model for given nominal pressure loss and nominal mass flow rate.
</p>
<p>
The pump controls a pressure ramp from 1.9 bar to 2.1 bar. 
This causes an appropriate ramp on the mass flow rate for pipe1, which has a boundary pressure of 1 bar. 
Flow reversal occurs in the pipe, which has a boundary pressure of 2 bar.
The Command plotResults can be used to see the pump speed N, which is controlled ideally to obtain the pressure ramp.
Moreover the Reynolds number as well as m_flow_turbulent and dp_turbulent are plotted for the pipe.
</p>
<p>
Next the pressureLoss model of the pipe could be investigated for <tt>length_nominal</tt>, 
which is obtained internally to fulfill the nominal pressure loss for given pipe diameter and roughness. 
Similarily the orifice could be investigated for <tt>zeta_nominal</tt>. 
Once the geometry has been designed, the NominalTurbulentFlow correlations can easily be replaced with 
WallFrictionPressureLoss correlations. Similarily the ControlledPump can be replaced with a PrescribedPump 
to investigate a real controller or with a Pump with rotational shaft to investigate inertia effects. 
</p>
</html>"));
end InverseParameterization;
