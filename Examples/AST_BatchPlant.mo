within Modelica_Fluid.Examples;
package AST_BatchPlant
  "Model of the experimental batch plant at Process Control Laboratory at University of Dortmund (Prof. Engell)"

  annotation (preferedView="info",Documentation(info="<html>
<p>
The process under consideration is an evaporation plant for a 
student lab at the Process Control Laboratory (AST) of the 
University of Dortmund that evaporates a water sodium chloride 
mixture so that a higher concentrated solution is produced. 
The task of the students is to learn how to program the process 
control system. A picture of the batch plant is shown in the figure
below.
</p>
 
<p align=\"center\">
<img src=\"../Images/Examples/AST_BatchPlant1.jpg\">
</p>
 
<p>
The flow sheet diagram is shown in the next figure.
</p>
 
<p align=\"center\">
<img src=\"../Images/Examples/AST_BatchPlant2.png\">
</p>
 
<p>
Pure water from tank B1 and concentrated sodium chloride 
solution from tank B2 are mixed in a mixing tank B3. 
After buffering in tank B4 the mixture flows to the 
evaporator B5. Here the water sodium chloride mixture 
is evaporated until the desired con-centration is reached. 
The steam is condensed in the condenser K1 and cooled 
afterwards in the cooling tank B6. The concentrated
 solution is also led to a cooling tank B7. The cooled 
fluids are pumped back to the charging vessels by the 
pumps P1 and P2. Be-tween the tanks several valves are 
present that are regulated by a central control system.
</p>
</html>"));

  model BatchPlant_StandardWater

    replaceable package BatchMedium = Modelica.Media.Water.StandardWater constrainedby
      Modelica.Media.Interfaces.PartialTwoPhaseMedium "Component media";

  /*
  replaceable package BatchMedium =Modelica.Media.Electrolytes.WaterNaCl extends 
    Modelica.Media.Interfaces.PartialTwoPhaseMedium "Component media";
*/

    parameter Modelica.SIunits.Length pipeDiameter = 0.01;

    BaseClasses.TankWith3InletOutletArraysWithEvaporatorCondensor B5(
      redeclare package Medium = BatchMedium,
      height=0.5,
      n_SidePorts=1,
      V0=0.001,
      n_BottomPorts=1,
      bottom_pipeArea={0.0001},
      top_pipeArea={0.0001},
      n_TopPorts=1,
      min_level_for_heating=0.0001,
      level_start=0.0009,
      crossArea=0.05,
      initType=BaseClasses.Init.InitialValues,
      side_pipeArea={0.0001}) 
      annotation (Placement(transformation(extent={{-100,-60},{-20,-20}},
            rotation=0)));
    Modelica_Fluid.Valves.ValveDiscrete V12(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={-80,2},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V15(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={-80,-82},
          extent={{10,-10},{-10,10}},
          rotation=90)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow HeatB5 
      annotation (Placement(transformation(extent={{-140,-50},{-120,-30}},
            rotation=0)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow CoolingB7 
      annotation (Placement(transformation(extent={{-140,-130},{-120,-110}},
            rotation=0)));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow CoolingB6 
      annotation (Placement(transformation(
          origin={110,-60},
          extent={{-10,-10},{10,10}},
          rotation=180)));

    Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.Controller controller(
        Transition3(enableTimer=true, waitTime=30), Transition7(
        condition=true,
        enableTimer=true,
        waitTime=300))       annotation (Placement(transformation(extent={{80,
              60},{120,100}}, rotation=0)));

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-300,
              -300},{300,300}}),      graphics),
      experiment(StopTime=3100),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/AST_BatchPlant_StandardWater/plot level.mos"
          "plot level"));
    Modelica_Fluid.Valves.ValveDiscrete V11(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(extent={{-60,78},{-40,98}}, rotation=
             0)));
    Modelica_Fluid.Valves.ValveDiscrete V8(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={-80,162},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V9(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={80,162},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V2(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(extent={{-30,230},{-50,250}},
            rotation=0)));
    Modelica_Fluid.Valves.ValveDiscrete V4(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(extent={{50,230},{30,250}}, rotation=
             0)));
    Modelica_Fluid.Valves.ValveDiscrete V3(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(extent={{-114,210},{-134,230}},
            rotation=0)));
    Fittings.TeeJunctionIdeal volume2(
      redeclare package Medium = BatchMedium) 
               annotation (Placement(transformation(
          origin={-180,220},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V6(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(extent={{112,210},{132,230}},
            rotation=0)));
    Fittings.TeeJunctionIdeal volume8(
      redeclare package Medium = BatchMedium) 
               annotation (Placement(transformation(
          origin={160,220},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V23(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={-110,-250},
          extent={{-10,10},{10,-10}},
          rotation=180)));
    Modelica_Fluid.Valves.ValveDiscrete V1(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={-180,110},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V22(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={-180,-56},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V5(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={160,110},
          extent={{10,-10},{-10,10}},
          rotation=270)));
    Modelica_Fluid.Valves.ValveDiscrete V24(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={90,-250},
          extent={{10,10},{-10,-10}},
          rotation=180)));
    Modelica_Fluid.Valves.ValveDiscrete V25(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={160,-20},
          extent={{10,-10},{-10,10}},
          rotation=270)));
    Modelica_Fluid.Valves.ValveDiscrete V20(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={60,-210},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V19(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={4,-200},
          extent={{-10,10},{10,-10}},
          rotation=180)));
    Modelica_Fluid.Valves.ValveDiscrete V10(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={20,70},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Modelica_Fluid.Valves.ValveDiscrete V21(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={30,-250},
          extent={{10,10},{-10,-10}},
          rotation=180)));
    Modelica_Fluid.Fittings.TeeJunctionVolume volume5(
      redeclare package Medium = BatchMedium,
      V=0.001) annotation (Placement(transformation(extent={{50,-260},{70,-240}},
            rotation=0)));
    Modelica_Fluid.Valves.ValveDiscrete V18(
      redeclare package Medium = BatchMedium,
      m_flow_nominal = 1,
      dp_nominal = 100) 
      annotation (Placement(transformation(
          origin={-56,-200},
          extent={{10,10},{-10,-10}},
          rotation=180)));
    Modelica_Fluid.Machines.PrescribedPumpNPSH P1(
      redeclare package Medium = BatchMedium,
      m_flow_start=0.1,
      redeclare function flowCharacteristic = 
          Modelica_Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow
          (                                                                   q_nominal={0,
              0.001,0.0015}, head_nominal={100,50,0}),
      N_nominal=200,
      use_N_input=true,
      V(displayUnit="ml") = 0.0001,
      checkValve=false,
      p_a_start=100000,
      p_b_start=100000) 
      annotation (Placement(transformation(extent={{-128,-260},{-148,-240}},
            rotation=0)));
    Modelica_Fluid.Machines.PrescribedPumpNPSH P2(
      redeclare package Medium = BatchMedium,
      m_flow_start=0.1,
      redeclare function flowCharacteristic = 
          Modelica_Fluid.Machines.BaseClasses.PumpCharacteristics.quadraticFlow
          (                                                                  q_nominal={0,
              0.001,0.0015}, head_nominal={100,50,0}),
      N_nominal=200,
      use_N_input=true,
      V(displayUnit="ml") = 0.0001,
      checkValve=true,
      p_a_start=100000,
      p_b_start=100000) 
      annotation (Placement(transformation(extent={{110,-260},{130,-240}},
            rotation=0)));
    Modelica_Fluid.Vessels.Tank B1(
      level_start=0.2,
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      crossArea=0.05,
      V0=0.0001,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)},
      stiffCharacteristicForEmptyPort=false) 
                         annotation (Placement(transformation(extent={{-100,180},
              {-60,220}}, rotation=0)));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{180,250},
              {200,270}},       rotation=0)));
    Modelica.Blocks.Logical.TriggeredTrapezoid P1_on(amplitude=100, rising=0) 
      annotation (Placement(transformation(extent={{-168,-230},{-148,-210}},
            rotation=0)));
    Modelica.Blocks.Logical.TriggeredTrapezoid P2_on(amplitude=50, rising=0) 
      annotation (Placement(transformation(extent={{90,-230},{110,-210}},
            rotation=0)));
    Modelica_Fluid.Vessels.Tank B2(
      level_start=0.2,
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      crossArea=0.05,
      V0=0.0001,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)},
      stiffCharacteristicForEmptyPort=false) 
                         annotation (Placement(transformation(extent={{60,180},
              {100,220}}, rotation=0)));
    Modelica_Fluid.Vessels.Tank B3(
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      crossArea=0.05,
      V0=0.0001,
      nTopPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0),Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=
           0.011, portLevel=0)},
      level_start=0.02,
      stiffCharacteristicForEmptyPort=false) 
                         annotation (Placement(transformation(extent={{-20,100},
              {20,140}}, rotation=0)));
    Vessels.Tank B4(
      redeclare package Medium = BatchMedium,
      levelMax=0.5,
      crossArea=0.05,
      V0=0.0001,
      level_start=0.015,
      nTopPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)},
      stiffCharacteristicForEmptyPort=false) 
                         annotation (Placement(transformation(extent={{-100,30},
              {-60,70}}, rotation=0)));
    Vessels.Tank B7(
      redeclare package Medium = BatchMedium,
      V0=0.0001,
      nTopPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)},
      level_start=0.009,
      levelMax=0.5,
      crossArea=0.05,
      stiffCharacteristicForEmptyPort=false,
      T_start=298,
      use_HeatTransfer=true,
      redeclare model HeatTransfer = 
          Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.ConstantHeatTransfer
          (alpha0=4.9e3)) 
                         annotation (Placement(transformation(extent={{-100,
              -140},{-60,-100}}, rotation=0)));
    Pipes.DistributedPipe pipeB1B2(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=pipeDiameter,
      height_ab=0) annotation (Placement(transformation(extent={{10,230},{-10,
              250}}, rotation=0)));
    Pipes.StaticPipe pipeB1B3(
      redeclare package Medium = BatchMedium,
      length=1,
      redeclare model PressureLoss = 
          Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedWallFriction,
      diameter=100*pipeDiameter,
      height_ab=0.1) 
      annotation (Placement(transformation(extent={{-42,134},{-62,154}},
            rotation=0)));
    Pipes.StaticPipe pipeB2B3(
      redeclare package Medium = BatchMedium,
      length=1,
      redeclare model PressureLoss = 
          Modelica_Fluid.Pipes.BaseClasses.PressureLoss.DetailedWallFriction,
      diameter=100*pipeDiameter,
      height_ab=0.1) 
      annotation (Placement(transformation(extent={{36,134},{56,154}}, rotation=
             0)));
    Pipes.DistributedPipe pipeB1B1(
      redeclare package Medium = BatchMedium,
      height_ab=0.5,
      diameter=pipeDiameter,
      length=1,
      nNodes=1,
      modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb) 
      annotation (Placement(transformation(
          origin={20,30},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Pipes.StaticPipe pipeB6Pump(
      redeclare package Medium = BatchMedium,
      length=0.5,
      diameter=pipeDiameter,
      height_ab=0.5) 
                   annotation (Placement(transformation(
          origin={60,-106},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Pipes.StaticPipe pipeB7Pump(
      redeclare package Medium = BatchMedium,
      length=1,
      diameter=pipeDiameter,
      height_ab=0.1) 
                   annotation (Placement(transformation(
          origin={-80,-170},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Pipes.DistributedPipe pipePump1B1(
      redeclare package Medium = BatchMedium,
      diameter=pipeDiameter,
      height_ab=3,
      length=3,
      nNodes=1,
      modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb) 
                   annotation (Placement(transformation(
          origin={-180,30},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Pipes.DistributedPipe pipePump2B2(
      redeclare package Medium = BatchMedium,
      diameter=pipeDiameter,
      height_ab=3,
      length=3,
      nNodes=1,
      modelStructure=Modelica_Fluid.Types.ModelStructure.a_vb) 
                   annotation (Placement(transformation(
          origin={160,30},
          extent={{-10,10},{10,-10}},
          rotation=90)));
    Vessels.Tank B6(
      redeclare package Medium = BatchMedium,
      V0=0.0001,
      nTopPorts=1,
      levelMax=0.5,
      crossArea=0.05,
      level_start=0.02,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.011,
          portLevel=0)},
      stiffCharacteristicForEmptyPort=false,
      T_start=298,
      use_HeatTransfer=true,
      redeclare model HeatTransfer = 
          Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.ConstantHeatTransfer
          (alpha0=4.9e3)) 
                         annotation (Placement(transformation(extent={{80,-80},
              {40,-40}}, rotation=0)));
    Fittings.MultiPort multiPort(redeclare package Medium = BatchMedium,
        nPorts_b=3) annotation (Placement(transformation(
          extent={{4,10},{-4,-10}},
          rotation=-90,
          origin={-26,-228})));
    Modelica_Fluid.Fittings.TeeJunctionVolume volume4(
      redeclare package Medium = BatchMedium,
      V=0.001) annotation (Placement(transformation(extent={{-36,-260},{-16,
              -240}},
            rotation=0)));
  equation
    controller.sensors.LIS_301 = B3.level;
    controller.sensors.QI_302 = 0;//B3.medium.X[2];
    controller.sensors.LIS_501 = B5.level;
    controller.sensors.QIS_502 = 0;//B5.medium.X[2];
    controller.sensors.TI_503 = B5.medium.T;
    controller.sensors.LIS_601 = B6.level;
    controller.sensors.TIS_602 = B6.medium.T;
    controller.sensors.LIS_701 = B7.level;
    controller.sensors.TIS_702 = B7.medium.T;

    P1_on.u = controller.actuators.P1;
    P2_on.u = controller.actuators.P2;
    V1.open = controller.actuators.V1;
    V2.open = controller.actuators.V2;
    V3.open = controller.actuators.V3;
    V4.open = controller.actuators.V4;
    V5.open = controller.actuators.V5;
    V6.open = controller.actuators.V6;
    V8.open = controller.actuators.V8;
    V9.open = controller.actuators.V9;
    V10.open = controller.actuators.V10;
    V11.open = controller.actuators.V11;
    V12.open = controller.actuators.V12;
    V15.open = controller.actuators.V15;
    V18.open = controller.actuators.V18;
    V19.open = controller.actuators.V19;
    V20.open = controller.actuators.V20;
    V21.open = controller.actuators.V21;
    V22.open = controller.actuators.V22;
    V23.open = controller.actuators.V23;
    V24.open = controller.actuators.V24;
    V25.open = controller.actuators.V25;
    HeatB5.Q_flow = if controller.actuators.T5_Heater then 20000 else 0;
    CoolingB6.Q_flow = if controller.actuators.T6_Cooling then -2000 else 0;
    CoolingB7.Q_flow = if controller.actuators.T7_Cooling then -2000 else 0;

    connect(V24.port_b, P2.port_a) annotation (Line(points={{100,-250},{110,
            -250}},                 color={0,127,255}));
    connect(V15.port_a, B5.BottomFluidPort[1]) annotation (Line(points={{-80,-72},
            {-80,-60.4}},      color={0,127,255}));
    connect(V3.port_a, B1.topPorts[1]) annotation (Line(points={{-114,220},{
            -106,220},{-106,230},{-79,230},{-79,220}}, color={0,127,255}));
    connect(B1.ports[1], V8.port_b) annotation (Line(points={{-81,180},{-81,176},
            {-80,176},{-80,172}},
          color={0,127,255}));
    connect(P1_on.y, P1.N_in) annotation (Line(points={{-147,-220},{-138,-220},
            {-138,-240}},     color={0,0,127}));
    connect(P2_on.y, P2.N_in) annotation (Line(points={{111,-220},{120,-220},{
            120,-240}},      color={0,0,127}));
    connect(B2.topPorts[1], V6.port_a) annotation (Line(points={{81,220},{81,
            228},{106,228},{106,220},{112,220}}, color={0,127,255}));
    connect(B2.ports[1], V9.port_b) annotation (Line(points={{79,180},{79,176},
            {80,176},{80,172}},
          color={0,127,255}));
    connect(V11.port_b, B3.ports[1]) annotation (Line(points={{-40,88},{-1,88},
            {-1,100}},
                    color={0,127,255}));
    connect(V10.port_b, B3.ports[2]) annotation (Line(points={{20,80},{20,88},{
            -1,88},{-1,100}},
                           color={0,127,255}));
    connect(B4.ports[1], V12.port_b) annotation (Line(points={{-81,30},{-81,21},
            {-80,21},{-80,12}},
          color={0,127,255}));
    connect(CoolingB7.port, B7.heatPort) annotation (Line(points={{-120,-120},{
            -100,-120}}, color={191,0,0}));
    connect(V2.port_a, pipeB1B2.port_b) annotation (Line(points={{-30,240},{-10,
            240}}, color={0,127,255}));
    connect(pipeB1B2.port_a, V4.port_b) annotation (Line(points={{10,240},{30,
            240}}, color={0,127,255}));
    connect(pipeB1B3.port_b, V8.port_a) annotation (Line(points={{-62,144},{-80,
            144},{-80,152}}, color={0,127,255}));
    connect(pipeB1B3.port_a, B3.topPorts[1]) annotation (Line(points={{-42,144},
            {1,144},{1,138}}, color={0,127,255}));
    connect(pipeB2B3.port_a, B3.topPorts[2]) annotation (Line(points={{36,144},
            {1,144},{1,142}}, color={0,127,255}));
    connect(pipeB2B3.port_b, V9.port_a) annotation (Line(points={{56,144},{80,
            144},{80,152}}, color={0,127,255}));
    connect(V11.port_a, B4.topPorts[1]) annotation (Line(points={{-60,88},{-79,
            88},{-79,70}}, color={0,127,255}));
    connect(V10.port_a, pipeB1B1.port_b) annotation (Line(points={{20,60},{20,
            40}}, color={0,127,255}));
    connect(B5.TopFluidPort[1], V12.port_a) annotation (Line(points={{-80,-19.6},
            {-80,-8}}, color={0,0,255}));
    connect(V15.port_b, B7.topPorts[1]) annotation (Line(points={{-80,-92},{-80,
            -96},{-80,-100},{-79,-100}},
                    color={0,127,255}));
    connect(B7.ports[1], pipeB7Pump.port_b) annotation (Line(points={{-81,-140},
            {-81,-150},{-80,-150},{-80,-160}},
                         color={0,127,255}));
    connect(pipePump1B1.port_a, V22.port_b) annotation (Line(points={{-180,20},
            {-180,-46}}, color={0,127,255}));
    connect(V25.port_b, pipePump2B2.port_a) annotation (Line(points={{160,-10},
            {160,20}},color={0,127,255}));
    connect(B6.ports[1], pipeB6Pump.port_b) annotation (Line(points={{61,-80},{
            61,-88},{60,-88},{60,-96}},
                      color={0,127,255}));
    connect(B6.topPorts[1], B5.Condensed) annotation (Line(points={{59,-40},{59,
            -28},{-19.6,-28}}, color={0,127,255}));
    connect(CoolingB6.port, B6.heatPort) annotation (Line(points={{100,-60},{80,
            -60}}, color={191,0,0}));
    connect(V20.port_b, pipeB6Pump.port_a) annotation (Line(points={{60,-200},{
            60,-116}}, color={0,127,255}));
    connect(HeatB5.port, B5.HeatPort) annotation (Line(points={{-120,-40},{-102,
            -40}}, color={191,0,0}));
    connect(volume2.port_3, V3.port_b)     annotation (Line(
        points={{-170,220},{-134,220}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V1.port_b, volume2.port_1)     annotation (Line(
        points={{-180,120},{-180,210}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume2.port_2, V2.port_b)     annotation (Line(
        points={{-180,230},{-180,240},{-50,240}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V6.port_b, volume8.port_3)     annotation (Line(
        points={{132,220},{150,220}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume8.port_2, V4.port_a)     annotation (Line(
        points={{160,230},{160,240},{50,240}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume8.port_1, V5.port_b)     annotation (Line(
        points={{160,210},{160,120}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume5.port_1, V21.port_b)     annotation (Line(
        points={{50,-250},{40,-250}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume5.port_2, V24.port_a)     annotation (Line(
        points={{70,-250},{80,-250}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume5.port_3, V20.port_a)     annotation (Line(
        points={{60,-240},{60,-220}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(P2.port_b, V25.port_a) annotation (Line(
        points={{130,-250},{160,-250},{160,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V22.port_a, P1.port_b) annotation (Line(
        points={{-180,-66},{-180,-250},{-148,-250}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V5.port_a, pipePump2B2.port_b) annotation (Line(
        points={{160,100},{160,40}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V1.port_a, pipePump1B1.port_b) annotation (Line(
        points={{-180,100},{-180,40}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(P1.port_a, V23.port_b) annotation (Line(
        points={{-128,-250},{-120,-250}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V23.port_a, volume4.port_1) annotation (Line(
        points={{-100,-250},{-36,-250}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume4.port_2, V21.port_a) annotation (Line(
        points={{-16,-250},{20,-250}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(multiPort.port_a, volume4.port_3) annotation (Line(
        points={{-26,-232},{-26,-240}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipeB1B1.port_a, multiPort.ports_b[2]) annotation (Line(
        points={{20,20},{20,-100},{-26,-100},{-26,-224}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(multiPort.ports_b[3], V19.port_b) annotation (Line(
        points={{-23.3333,-224},{-23.3333,-200},{-6,-200}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipeB7Pump.port_a, V18.port_a) annotation (Line(
        points={{-80,-180},{-80,-200},{-66,-200}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V18.port_b, multiPort.ports_b[1]) annotation (Line(
        points={{-46,-200},{-28.6667,-200},{-28.6667,-224}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(V19.port_a, pipeB6Pump.port_a) annotation (Line(
        points={{14,-200},{30,-200},{30,-160},{60,-160},{60,-116}},
        color={0,127,255},
        smooth=Smooth.None));
  end BatchPlant_StandardWater;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    block TriggeredTrapezoid "Triggered trapezoid generator"
      extends Modelica.Blocks.Interfaces.partialBooleanBlockIcon;

      parameter Real amplitude=1 "Amplitude of trapezoid";
      parameter Modelica.SIunits.Time rising(final min=0)=0
        "Rising duration of trapezoid";
      parameter Modelica.SIunits.Time falling(final min=0)=rising
        "Falling duration of trapezoid";
      parameter Real offset=0 "Offset of output signal";

      Modelica.Blocks.Interfaces.BooleanInput u
        "Connector of Boolean input signal" 
                                       annotation (Placement(transformation(
              extent={{-140,-20},{-100,20}}, rotation=0)));
      Modelica.Blocks.Interfaces.RealOutput y "Connector of Real output signal"
        annotation (Placement(transformation(extent={{100,-10},{120,10}},
              rotation=0)));

      annotation (
        Icon(graphics={
            Line(points={{-60,-70},{-60,-70},{-30,40},{8,40},{40,-70},{40,-70}}),
            Line(points={{-90,-70},{82,-70}}, color={192,192,192}),
            Line(points={{-80,68},{-80,-80}}, color={192,192,192}),
            Polygon(
              points={{90,-70},{68,-62},{68,-78},{90,-70}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-80,90},{-88,68},{-72,68},{-80,90}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,-70},{-60,-70},{-60,24},{8,24},{8,-70},{60,-70}},
                color={255,0,255})}),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Line(points={{-80,-20},{-60,-20},{-30,40},{8,40},{40,-20},{60,-20}}),
            Line(points={{-90,-70},{82,-70}}, color={0,0,0}),
            Line(points={{-80,68},{-80,-80}}, color={0,0,0}),
            Polygon(
              points={{90,-70},{68,-62},{68,-78},{90,-70}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-80,90},{-88,68},{-72,68},{-80,90}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,-68},{-60,-68},{-60,-42},{8,-42},{8,-68},{60,-68}},
                color={255,0,255}),
            Line(
              points={{-60,40},{-60,-42}},
              color={0,0,0},
              pattern=LinePattern.Dot),
            Line(
              points={{8,-42},{8,40}},
              color={0,0,0},
              pattern=LinePattern.Dot),
            Line(points={{-58,40},{-28,40}}, color={0,0,0}),
            Line(points={{8,-20},{40,-20}}, color={0,0,0}),
            Line(points={{-20,40},{-20,-20}}, color={0,0,0}),
            Line(points={{-20,-20},{-20,-70}}, color={0,0,0}),
            Text(
              extent={{-42,48},{-42,38}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              textString="rising"),
            Text(
              extent={{24,-10},{24,-20}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              textString="falling"),
            Polygon(
              points={{-58,40},{-54,42},{-54,38},{-58,40}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-30,40},{-34,42},{-34,38},{-30,40}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{8,-20},{12,-18},{12,-22},{8,-20}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{40,-20},{36,-18},{36,-22},{40,-20}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-22,-24},{-20,-20},{-18,-24},{-22,-24}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-18,-66},{-22,-66},{-20,-70},{-18,-66}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-22,36},{-20,40},{-18,36},{-22,36}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{-18,-16},{-22,-16},{-20,-20},{-18,-16}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-40,6},{0,-4}},
              lineColor={255,255,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-20,6},{-20,-4}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              textString="amplitude"),
            Rectangle(
              extent={{-40,-48},{0,-58}},
              lineColor={255,255,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-20,-48},{-20,-58}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              textString="offset"),
            Text(
              extent={{60,-82},{94,-92}},
              lineColor={0,0,0},
              textString="time"),
            Text(
              extent={{-88,-4},{-54,-14}},
              lineColor={0,0,0},
              textString="y"),
            Text(
              extent={{-88,-46},{-54,-56}},
              lineColor={0,0,0},
              textString="u"),
            Polygon(
              points={{40,60},{36,62},{36,58},{40,60}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-56,60},{36,60}}, color={0,0,0}),
            Polygon(
              points={{-60,60},{-56,62},{-56,58},{-60,60}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-16,72},{-16,62}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid,
              textString="y_high")}),
        Documentation(info="<HTML>
<p>The block TriggeredTrapezoid has a boolean input and a real
output signal and requires the parameters <i>amplitude</i>,
<i>rising</i>, <i>falling</i> and <i>offset</i>. The
output signal <b>y</b> represents a trapezoidal signal dependent on the
input signal <b>u</b>.
</p>
<p>The behaviour is as follows: Assume the initial input to be false. In this
case, the output will be <i>offset</i>. After a rising edge (i.e. the input
changes from false to true), the output is rising during <i>rising</i> to the
sum of <i>offset</i> and <i>amplitude</i>. In contrast, after a falling
edge (i.e. the input changes from true to false), the output is falling
during <i>falling</i> to a value of <i>offset</i>.
</p>
<p>Note, that the case of edges before expiration of rising or falling is
handled properly.</p>
</HTML>
"));

    protected
      discrete Real endValue "Value of y at time of recent edge";
      discrete Real rate "Current rising/falling rate";
      discrete Modelica.SIunits.Time T
        "Predicted time of output reaching endValue";
    public
      Modelica.Blocks.Interfaces.BooleanOutput y_high 
        annotation (Placement(transformation(extent={{100,-90},{120,-70}},
              rotation=0)));
    initial equation
      /* A start value of y is set, because pre(y) is present
     to avoid a warning message from the compiler. However,
     this setting does not have an effect, because y is initialized
     correctly, before pre(y) is used
  */
      pre(y) = 0;
    equation
        y_high = time < T;
        y = if y_high then endValue - (T - time)*rate else  endValue;

        when {initial(),u,not u} then
          endValue = if u then offset + amplitude else offset;
          rate = if u and (rising > 0) then amplitude/rising else 
            if not u and (falling > 0) then -amplitude/falling else 0;
          T = if u and not (rising > 0) or not u and not (falling
             > 0) or not abs(amplitude) > 0 or initial() then time else time
             + (endValue - pre(y))/rate;
        end when;
    end TriggeredTrapezoid;

    block setReal "Set output signal to a time varying Real expression"

      Modelica.Blocks.Interfaces.RealInput u "Set value of Real input" 
        annotation (                            Dialog(group=
              "Time varying input signal"), Placement(transformation(extent={{
                -140,-20},{-100,20}}, rotation=0)));

      annotation (
        Window(
          x=0.29,
          y=0.23,
          width=0.6,
          height=0.6),
        Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={
            Rectangle(
              extent={{-100,40},{100,-40}},
              lineColor={0,0,0},
              fillColor={235,235,235},
              fillPattern=FillPattern.Solid,
              borderPattern=BorderPattern.Raised),
            Text(
              extent={{-96,15},{96,-15}},
              lineColor={0,0,0},
              fillColor={0,255,0},
              fillPattern=FillPattern.Solid,
              textString="%u"),
            Text(extent={{-150,90},{140,50}}, textString="%name")}),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics),
        Documentation(info="<html>
 
</html>"));

    end setReal;

    model ValveDiscrete "Valve for water/steam flows with linear pressure drop"
      extends Modelica_Fluid.Fittings.BaseClasses.PartialTwoPortTransport;
      parameter Modelica_Fluid.Types.HydraulicConductance Kv
        "Hydraulic conductance at full opening";
      Modelica.Blocks.Interfaces.BooleanInput open annotation (Placement(
            transformation(
            origin={0,84},
            extent={{10,-10},{-10,10}},
            rotation=90)));
      parameter Real m_flow_small = 1e-7 "massflow when valve is closed";
      parameter Boolean finiteRiseTime = false
        "= true, if valve is opened/closed linearly in riseTime time";
      parameter Real riseTime = 0 "Time to open or close the valve" 
                                           annotation(Dialog(enable=finiteRiseTime));
    protected
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.TriggeredTrapezoid
        trapezoid(rising=riseTime) if finiteRiseTime 
        annotation (Placement(transformation(
            origin={0,50},
            extent={{-10,-10},{10,10}},
            rotation=270)));
      Modelica.Blocks.Sources.RealExpression m_flow_trapezoid(y=if trapezoid.y_high then 
                  trapezoid.y*Kv*dp else Kv*dp*m_flow_small) if finiteRiseTime 
        annotation (Placement(transformation(extent={{-90,-40},{40,-20}},
              rotation=0)));
      Modelica.Blocks.Sources.RealExpression m_flow_pulse(y=if open then Kv*dp else 
                  Kv*dp*m_flow_small) if not finiteRiseTime 
        annotation (Placement(transformation(extent={{-90,-60},{40,-40}},
              rotation=0)));
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.setReal set(u=m_flow) 
        annotation (Placement(transformation(extent={{70,-50},{90,-30}},
              rotation=0)));
    annotation (
      Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Line(points={{0,50},{0,0}}, color={0,0,0}),
            Rectangle(
              extent={{-20,60},{20,50}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Text(extent={{-145,-58},{146,-98}}, textString="%name"),
            Polygon(
              points={{-100,50},{100,-50},{100,50},{0,0},{-100,-50},{-100,50}},
              fillColor=DynamicSelect({255,255,255}, if open > 0.5 then {0,255,
                  0} else {255,255,255}),
              lineColor={0,0,0})}),
      Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics),
      Documentation(info="<HTML>
<p>This very simple model provides a pressure drop which is proportional to the flowrate and to the <tt>opening</tt> signal, without computing any fluid property.
<p>A medium model must be nevertheless be specified, so that the fluid ports can be connected to other components using the same medium model.
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));

    equation
      connect(trapezoid.u, open) annotation (Line(points={{2.20436e-015,62},{0,
              62},{0,84}}, color={255,0,255}));
      connect(m_flow_trapezoid.y, set.u) annotation (Line(points={{46.5,-30},{
              60,-30},{60,-40},{68,-40}}, color={0,0,127}));
      connect(m_flow_pulse.y, set.u) annotation (Line(points={{46.5,-50},{60,
              -50},{60,-40},{68,-40}}, color={0,0,127}));
    end ValveDiscrete;

    model TankWith3InletOutletArraysWithEvaporatorCondensor
      "Tank with Heating and Evaporation"
      import Modelica.SIunits.Conversions.*;
      import Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.Init;
      replaceable package Medium = Modelica.Media.Water.WaterIF97_ph 
        constrainedby Modelica.Media.Interfaces.PartialTwoPhaseMedium
        "Medium in the component" 
        annotation (choicesAllMatching=true);
    // parameter for Tank
      parameter Modelica.SIunits.Area crossArea "Tank area";
      parameter SI.Area top_pipeArea[n_TopPorts] "Area of outlet pipe";
      parameter SI.Area side_pipeArea[n_SidePorts] "Area of outlet pipe";
      parameter SI.Area bottom_pipeArea[n_BottomPorts] "Area of outlet pipe";
      parameter Modelica.SIunits.Height height(min=0) = 10 "Height of Tank";
      parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";
      constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
      parameter Real side_heights[n_SidePorts]=zeros(n_SidePorts);
      parameter Real bottom_heights[n_BottomPorts]=zeros(n_BottomPorts);
      parameter Real top_heights[n_TopPorts]=fill(height, n_TopPorts);
      parameter SI.Height level_start(min=0) "Initial tank level" 
        annotation(Dialog(tab="Initialization"));
      parameter Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.Init
        initType =                                 Init.GuessValues
        "Initialization option" 
        annotation(Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true
        "Use T_start if true, otherwise h_start"                                    annotation(Dialog(tab = "Initialization"), Evaluate = true);
      parameter Medium.Temperature T_start=if use_T_start then 293.15 else 
          Medium.temperature_phX(p_ambient, h_start, X_start)
        "Start value of temperature" 
        annotation(Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=if use_T_start then Medium.specificEnthalpy_pTX(
          p_ambient, T_start, X_start[1:Medium.nXi]) else 1e4
        "Start value of specific enthalpy" 
        annotation(Dialog(tab = "Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX]=Medium.reference_X
        "Start value of mass fractions m_i/m" 
        annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
      parameter Medium.AbsolutePressure p_ambient=101325
        "Tank surface pressure";
      parameter Medium.Temperature T_ambient=293.15 "Tank surface Temperature";
      parameter Integer n_TopPorts=1 "number of Top connectors";
      parameter Integer n_SidePorts=1 "number of side connectors";
      parameter Integer n_BottomPorts=1 "number of bootom connectors";
      Medium.BaseProperties medium(
        preferredMediumStates=true,
        p(start=p_ambient),
        T(start=T_start),
        Xi(start=X_start[1:Medium.nXi]));
      Modelica.SIunits.Height level(
        stateSelect=StateSelect.prefer,
        min=0,
        max=height) "Level height of tank";
      SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
      SI.Energy U "Internal energy of tank volume";
      Real m(quantity=Medium.mediumName, unit="kg") "Mass of tank volume";
      Real mXi[Medium.nXi](quantity=Medium.substanceNames, each unit="kg")
        "Component masses of the independent substances";
    // Hilfsvariablen
      Real H_flow_BottomPorts[n_BottomPorts];
      Real H_flow_SidePorts[n_SidePorts];
      Real H_flow_TopPorts[n_TopPorts];
      Real m_flow_BottomPorts[n_BottomPorts];
      Real m_flow_SidePorts[n_SidePorts];
      Real m_flow_TopPorts[n_TopPorts];

      Real m_flow_BottomPorts_pos[n_BottomPorts];
      Real m_flow_SidePorts_pos[n_SidePorts];
      Real m_flow_TopPorts_pos[n_TopPorts];
      Real m_flow_pos;
      Medium.MassFlowRate mXi_flow_topPorts[n_TopPorts,Medium.nXi];
      Medium.MassFlowRate mXi_flowBottomPorts[n_BottomPorts,Medium.nXi];
      Medium.MassFlowRate mXi_flow_sidePorts[n_SidePorts,Medium.nXi];

    // Connectors and InnerTanks
      Modelica_Fluid.Interfaces.FluidPort_b BottomFluidPort[n_BottomPorts](
        redeclare package Medium = Medium,
        m_flow(each start=0)) 
        annotation (Placement(transformation(
            origin={-100,-102},
            extent={{-10,-10},{10,10}},
            rotation=90)));
      Modelica_Fluid.Interfaces.FluidPort_a TopFluidPort[n_TopPorts](
        redeclare package Medium = Medium,
        m_flow(each start=0)) 
        annotation (Placement(transformation(extent={{-110,92},{-90,112}},
              rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_b SideFluidPort[n_SidePorts](
        redeclare package Medium = Medium,
        m_flow(each start=0)) 
        annotation (Placement(transformation(extent={{0,-10},{20,10}}, rotation=
               0)));
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.InnerTank InnerTankTop[n_TopPorts](
        each h=medium.h,
        each p_ambient=p_ambient,
        each d=medium.d,
        each Xi = medium.Xi,
        aboveLevel={level - top_heights[i] for i in 1:n_TopPorts},
        pipeArea={top_pipeArea[i] for i in 1:n_TopPorts},
        redeclare package Medium = Medium) 
          annotation (Placement(transformation(extent={{-140,60},{-120,80}},
              rotation=0)));
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.InnerTank
        InnerTankSide[                                                           n_SidePorts](
        each h=medium.h,
        each p_ambient=p_ambient,
        each d=medium.d,
        each Xi = medium.Xi,
        aboveLevel={level - side_heights[i] for i in 1:n_SidePorts},
        pipeArea={side_pipeArea[i] for i in 1:n_SidePorts},
        redeclare package Medium = Medium) 
          annotation (Placement(transformation(extent={{-20,0},{0,20}},
              rotation=0)));
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.InnerTank
        InnerTankBottom[                                                           n_BottomPorts](
        each h=medium.h,
        each p_ambient=p_ambient,
        each d=medium.d,
        each Xi = medium.Xi,
        aboveLevel={level - bottom_heights[i] for i in 1:n_BottomPorts},
        pipeArea={bottom_pipeArea[i] for i in 1:n_BottomPorts},
        redeclare package Medium = Medium) 
          annotation (Placement(transformation(extent={{-80,-80},{-60,-60}},
              rotation=0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPort 
        annotation (Placement(transformation(extent={{-220,-10},{-200,10}},
              rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_b Condensed(redeclare package Medium
          =        Medium) 
        annotation (Placement(transformation(extent={{192,50},{212,70}},
              rotation=0)));
    // parameter for Evaporator
      parameter Real min_level_for_heating;
      parameter Real k=4.9 "Wärmedurchgangskoeffizient";
      Real Q_lost "Wärmeverlust";

      Medium.SaturationProperties sat
        "State vector to compute saturation properties";
      Medium.SpecificEnthalpy h_v=Medium.dewEnthalpy(sat)
        "specific enthalpy of vapour";
      Medium.SpecificEnthalpy h_l=Medium.bubbleEnthalpy(sat)
        "specific enthalpy of liquid";
      Medium.SpecificEnthalpy h "'is'specific enthalpy of liquid";
      Medium.Density rho_v=Medium.dewDensity(sat) "density in vapour phase";
      Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
      Medium.Density rho "'is' density in liquid phase";

    equation
      H_flow_TopPorts   = InnerTankTop.H_flow;
      m_flow_TopPorts   = InnerTankTop.m_flow;
      mXi_flow_topPorts = InnerTankTop.mXi_flow;

      H_flow_SidePorts   = InnerTankSide.H_flow;
      m_flow_SidePorts   = InnerTankSide.m_flow;
      mXi_flow_sidePorts = InnerTankSide.mXi_flow;

      H_flow_BottomPorts  = InnerTankBottom.H_flow;
      m_flow_BottomPorts  = InnerTankBottom.m_flow;
      mXi_flowBottomPorts = InnerTankBottom.mXi_flow;

      for i in 1:n_BottomPorts loop
        m_flow_BottomPorts_pos[i] = (if m_flow_BottomPorts[i] > 0 then 
          m_flow_BottomPorts[i] else 0);
      end for;
      for i in 1:n_SidePorts loop
        m_flow_SidePorts_pos[i] = if m_flow_SidePorts[i] > 0 then m_flow_SidePorts[
          i] else 0;
      end for;
      for i in 1:n_TopPorts loop
        m_flow_TopPorts_pos[i] = if m_flow_TopPorts[i] > 0 then m_flow_TopPorts[i] else 
                0;
      end for;
      for i in 1:n_BottomPorts loop
        connect(InnerTankBottom[i].port, BottomFluidPort[i]) annotation (Line(
              points={{-70,-81},{-70,-102},{-100,-102}}, color={0,0,255}));
      end for;
      for i in 1:n_TopPorts loop
        connect(InnerTankTop[i].port, TopFluidPort[i])  annotation (Line(points=
               {{-130,59},{-92,59},{-92,102},{-100,102}}, color={0,0,255}));
      end for;
      for i in 1:n_SidePorts loop
        connect(InnerTankSide[i].port, SideFluidPort[i])  annotation (Line(
              points={{-10,-1},{24,-1},{24,0},{10,0}}, color={0,0,255}));
      end for;

      medium.p = p_ambient;
      medium.T = HeatPort.T;
    // Mass balance
      der(m) = sum(m_flow_BottomPorts) + sum(m_flow_SidePorts) + sum(
        m_flow_TopPorts) + Condensed.m_flow;
    // Energy balance

      U = m*medium.h - p_ambient*V "Internal energy of fluid";
      Q_lost = -k*(2*crossArea + 2*sqrt(Modelica.Constants.pi*crossArea))*level*(medium.T -
        T_ambient);
      m = V*medium.d "Mass of fluid";
      V = crossArea*level + V0 "Volume of fluid";
      mXi = m*medium.Xi "Mass of fluid components";
      sat.psat = medium.p;
      sat.Tsat = Medium.saturationTemperature(medium.p);

      if noEvent(medium.T < sat.Tsat) then
        if Medium.singleState then
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(
            H_flow_TopPorts) + Condensed.m_flow*actualStream(Condensed.h_outflow) + HeatPort.Q_flow + Q_lost
            "Mechanical work is neglected";
        else
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(
            H_flow_TopPorts) + Condensed.m_flow*actualStream(Condensed.h_outflow) - p_ambient*der(V) + Q_lost +
            HeatPort.Q_flow;
        end if;
        Condensed.h_outflow = h;
        Condensed.m_flow = 0;
        rho = medium.d;
        h = medium.h;
      else
        if Medium.singleState then
                                 //Q_flow cooling = - (HeatPort.Q_flow-Q_lost)
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(H_flow_TopPorts)
             + Condensed.m_flow*actualStream(Condensed.h_outflow)
            "Mechanical work is neglected";
        else
          der(U) = sum(H_flow_BottomPorts) + sum(H_flow_SidePorts) + sum(
            H_flow_TopPorts) + Condensed.m_flow*actualStream(Condensed.h_outflow) - p_ambient*der(V);
        end if;
        Condensed.h_outflow = h;
        Condensed.m_flow = -(HeatPort.Q_flow - Q_lost)/(h_v - h_l);
        rho = rho_l;//Density = liquid Densety
        h = h_l;    //Enthalpy = liquid Enthalpy
        if noEvent(HeatPort.Q_flow > 0.0) then
          assert(noEvent(abs(m_flow_pos) <= 0.01), "Es wird beim verdampfen befüllt.");
        end if;

      end if;

      m_flow_pos = sum(m_flow_TopPorts_pos) + sum(m_flow_SidePorts_pos) + sum(
        m_flow_BottomPorts_pos);

      for i in 1:Medium.nXi loop
           der(mXi[i]) = sum(mXi_flowBottomPorts[:,i]) +
                         sum(mXi_flow_sidePorts[:,i]) +
                         sum(mXi_flow_topPorts[:,i]);
      end for;

      assert(level < height, " 
    Tank ist überfüllt.
    ");

      assert(not (HeatPort.Q_flow > 0.0 and level <= min_level_for_heating), "
    Es wird leerer Tank bezeizt
  ");

    initial equation
      if initType == Init.GuessValues then
        // no initial equations
      elseif initType == Init.InitialValues then
        level = level_start;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      elseif initType == Init.SteadyStateHydraulic then
        der(level) = 0;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      else
        assert(false, "Unsupported initialization option");
      end if;
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},{
                200,100}}), graphics={
            Rectangle(
              extent={{-200,100},{0,-90}},
              lineColor={255,255,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent=DynamicSelect({{-200,-100},{0,0}}, {{-200,-100},{0,(-100
                   + 200*level/levelMax)}}),
              lineColor={0,127,255},
              fillColor={85,170,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-200,100},{-200,-100},{0,-100},{0,100}}, color={0,0,0}),
            Text(
              extent={{-198,74},{0,38}},
              lineColor={0,0,255},
              fillColor={0,127,255},
              fillPattern=FillPattern.Solid,
              textString="%name"),
            Text(
              extent={{-184,-64},{-14,-86}},
              lineColor={0,0,0},
              textString="%level_start"),
            Text(
              extent={{-192,-34},{-12,-54}},
              lineColor={0,0,0},
              textString="level_start ="),
            Line(
              points={{-200,100},{0,100}},
              color={0,0,0},
              pattern=LinePattern.Dot),
            Polygon(
              points={{0,100},{200,70},{200,50},{200,50},{0,80},{0,100}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,0,255}),
            Polygon(
              points={{20,98},{30,74},{52,84},{66,72},{86,78},{98,66},{118,74},
                  {130,60},{144,70},{152,60},{168,66},{180,54},{196,74},{190,76},
                  {180,64},{170,70},{156,66},{148,76},{132,68},{120,80},{100,74},
                  {88,88},{70,78},{50,92},{32,82},{28,100},{20,98},{20,98}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={170,255,255})}),
        Documentation(info="<HTML>
<p>
<p>This tank has the same geometric variables as TankWith3InletOutletArrays plus the feature of a HeatPort and the possibility of evaporation. 
(Assumption: The gas is condensed emidiatly afterwards so that a liquid boiling fluid is created.)
<p>The tank can be initialized with the following options:
<ul>
<li>GuessValues: no explicit initial conditions
<li>InitialValues: initial values of temperature (or specific enthalpy), composition and level are specified
<li>SteadyStateHydraulic: initial values of temperature (or specific enthalpy) and composition are specified; the initial level is determined so that levels and pressure are at steady state.
</ul>
Full steady state initialization is not supported, because the corresponding intial equations for temperature/enthalpy are undetermined (the flow rate through the port at steady state is zero). 
</p>
</HTML>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,-100},
                {200,100}}), graphics));

    end TankWith3InletOutletArraysWithEvaporatorCondensor;

    model InnerTank
        replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching=true);

        Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = 
            Medium) 
        annotation (Placement(transformation(
            origin={0,-110},
            extent={{-10,-10},{10,10}},
            rotation=90)));
        Boolean m_flow_negative( start = true) "true= massflow out of tank";
        constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
        input Real aboveLevel;
        input Real d;
        input Real p_ambient;
        input Real h;
        input Medium.MassFraction Xi[Medium.nXi]
        "Actual mass fractions of fluid in tank"                    annotation(Dialog);
        input Real pipeArea;
        output Real H_flow;
        output Real m_flow;
       output Medium.MassFlowRate mXi_flow[Medium.nXi]
        "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";

    equation
    m_flow_negative = (pre(m_flow_negative) and not port.p>p_ambient) or (port.m_flow < -1e-6);

    if noEvent(aboveLevel > 0) then
      port.p = aboveLevel*g*d + p_ambient - smooth(2,noEvent(if noEvent(m_flow < 0) then m_flow^2/(2*d*pipeArea^2) else 0));
    else
     if pre(m_flow_negative) then
        port.m_flow = 0;
      else
        port.p = p_ambient;
      end if;
    end if;

      H_flow = port.m_flow*actualStream(port.h_outflow);
      m_flow = port.m_flow;
      mXi_flow = port.m_flow*actualStream(port.Xi_outflow);
      port.h_outflow = h;
      port.Xi_outflow = Xi;
    end InnerTank;

    model Controller

      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ControllerUtilities.Port_Sensors
        sensors 
        annotation (Placement(transformation(extent={{-280,-40},{-200,40}},
              rotation=0)));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-200,
                -200},{200,200}}),      graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-200,-200},{
                200,200}}), graphics={
            Rectangle(
              extent={{-200,200},{200,-200}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-288,286},{262,208}},
              lineColor={0,0,255},
              textString="%name"),
            Line(points={{-48,0},{0,0}}, color={0,0,0}),
            Rectangle(extent={{-170,60},{-50,-60}}, lineColor={0,0,0}),
            Line(points={{0,40},{0,-40}}, color={0,0,0}),
            Line(points={{0,0},{26,0}}, color={0,0,0}),
            Polygon(
              points={{26,10},{50,0},{26,-10},{26,10}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Rectangle(extent={{50,60},{170,-60}}, lineColor={0,0,0}),
            Polygon(
              points={{-24,10},{0,0},{-24,-10},{-24,10}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid)}));
      Modelica_Fluid.Examples.AST_BatchPlant.BaseClasses.ControllerUtilities.Port_Actuators
        actuators 
        annotation (Placement(transformation(extent={{200,-20},{240,20}},
              rotation=0)));

      parameter Real w_dilution=0.003;
      parameter Real w_concentrate=0.005;
      parameter Real startTime=1;
      parameter Real T5_batch_level=0.211;

      Modelica.StateGraph.InitialStep InitialStep1 
        annotation (Placement(transformation(extent={{-180,90},{-160,110}},
              rotation=0)));
      Modelica.StateGraph.Transition Transition1(enableTimer=true, waitTime=
            startTime) annotation (Placement(transformation(extent={{-150,90},{
                -130,110}}, rotation=0)));
      Modelica.StateGraph.Step Step1 
        annotation (Placement(transformation(extent={{-120,90},{-100,110}},
              rotation=0)));
      Modelica.StateGraph.Transition Transition2(condition=LIS_301 >= 0.13) 
        annotation (Placement(transformation(extent={{-90,90},{-70,110}},
              rotation=0)));
      Modelica.StateGraph.Step Step2 annotation (Placement(transformation(
              extent={{-60,90},{-40,110}}, rotation=0)));
      Modelica.StateGraph.Transition Transition3(
        condition=true,
        enableTimer=true,
        waitTime=500) 
        annotation (Placement(transformation(extent={{-30,90},{-10,110}},
              rotation=0)));
      Modelica.StateGraph.Step Step3 annotation (Placement(transformation(
              extent={{0,90},{20,110}}, rotation=0)));
      Modelica.StateGraph.Transition Transition4(condition=LIS_301 <= 0.01) 
        annotation (Placement(transformation(extent={{30,90},{50,110}},
              rotation=0)));
      Modelica.StateGraph.Step Step4 annotation (Placement(transformation(
              extent={{60,90},{80,110}}, rotation=0)));
      Modelica.StateGraph.Transition Transition5(condition=T5_idle) 
        annotation (Placement(transformation(extent={{90,90},{110,110}},
              rotation=0)));
      Modelica.StateGraph.Step Step5 annotation (Placement(transformation(
              extent={{120,90},{140,110}}, rotation=0)));
      Modelica.StateGraph.Transition Transition6(condition=LIS_501 >=
            T5_batch_level) annotation (Placement(transformation(extent={{150,
                90},{170,110}}, rotation=0)));
      Modelica.StateGraph.Step Step6 annotation (Placement(transformation(
              extent={{-120,30},{-100,50}}, rotation=0)));
      Modelica.StateGraph.Transition Transition7(
        condition=true,
        enableTimer=true,
        waitTime=300) 
        annotation (Placement(transformation(extent={{-90,30},{-70,50}},
              rotation=0)));
      Modelica.StateGraph.Parallel Parallel1 annotation (Placement(
            transformation(extent={{-176,-100},{194,0}}, rotation=0)));
      Modelica.StateGraph.Step Step7 annotation (Placement(transformation(
              extent={{-122,-80},{-102,-60}}, rotation=0)));
      Modelica.StateGraph.Step Step8 annotation (Placement(transformation(
              extent={{-62,-80},{-42,-60}}, rotation=0)));
      Modelica.StateGraph.Step Step9 annotation (Placement(transformation(
              extent={{-2,-80},{18,-60}}, rotation=0)));
      Modelica.StateGraph.Step Step10 annotation (Placement(transformation(
              extent={{58,-80},{78,-60}}, rotation=0)));
      Modelica.StateGraph.Step Step11 annotation (Placement(transformation(
              extent={{118,-80},{138,-60}}, rotation=0)));
      Modelica.StateGraph.Step Step12 annotation (Placement(transformation(
              extent={{-62,-40},{-42,-20}}, rotation=0)));
      Modelica.StateGraph.Step Step13 annotation (Placement(transformation(
              extent={{-2,-40},{18,-20}}, rotation=0)));
      Modelica.StateGraph.Step Step14 annotation (Placement(transformation(
              extent={{58,-40},{78,-20}}, rotation=0)));
      Modelica.StateGraph.Transition Transition8(condition=T7_idle) 
        annotation (Placement(transformation(extent={{-92,-80},{-72,-60}},
              rotation=0)));
      Modelica.StateGraph.Transition Transition9(condition=LIS_501 <= 0.01) 
        annotation (Placement(transformation(extent={{-32,-80},{-12,-60}},
              rotation=0)));
      Modelica.StateGraph.Transition Transition10(condition=TIS_702 <= 298) 
        annotation (Placement(transformation(extent={{28,-80},{48,-60}},
              rotation=0)));
      Modelica.StateGraph.Transition Transition11(condition=LIS_701 <= 0.01) 
        annotation (Placement(transformation(extent={{88,-80},{108,-60}},
              rotation=0)));
      Modelica.StateGraph.Transition Transition12(condition=TIS_602 <= 298) 
        annotation (Placement(transformation(extent={{-32,-40},{-12,-20}},
              rotation=0)));
      Modelica.StateGraph.Transition Transition13(condition=LIS_601 <= 0.01) 
        annotation (Placement(transformation(extent={{28,-40},{48,-20}},
              rotation=0)));

      Real LIS_301;
      Real LIS_501;
      Real LIS_601;
      Real LIS_701;
      Real QI_302;
      Real QIS_502;
      Real TIS_602;
      Real TIS_702;
      Boolean T5_idle;
      Boolean T7_idle;
      Modelica.StateGraph.TransitionWithSignal TransitionWithSignal1 
        annotation (Placement(transformation(
            origin={-2,-150},
            extent={{-10,-10},{10,10}},
            rotation=180)));
      Modelica.Blocks.Sources.BooleanExpression BooleanExpression1(y=time >
            2500) 
        annotation (Placement(transformation(extent={{-104,-148},{-18,-116}},
              rotation=0)));
    equation
      LIS_301 = sensors.LIS_301;
      LIS_501 = sensors.LIS_501;
      LIS_601 = sensors.LIS_601;
      LIS_701 = sensors.LIS_701;
      QI_302 = sensors.QI_302;
      QIS_502 = sensors.QIS_502;
      TIS_602 = sensors.TIS_602;
      TIS_702 = sensors.TIS_702;
      T5_idle = not actuators.V12 and not actuators.V15 and not actuators.T5_Heater
         and sensors.LIS_501 < 0.01;
      T7_idle = not actuators.V15 and not actuators.V18 and not actuators.
        T7_Cooling and sensors.LIS_701 < 0.01;

      actuators.P1 = Step10.active;
      actuators.P2 = Step13.active;
      actuators.T5_Heater = Step6.active;
      actuators.T7_Cooling = Step9.active;
      actuators.T6_Cooling = Step12.active;
      actuators.V1 = Step10.active;
      actuators.V2 = false;
      actuators.V3 = Step10.active;
      actuators.V4 = false;
      actuators.V5 = Step13.active;
      actuators.V6 = Step13.active;
      actuators.V8 = Step1.active;
      actuators.V9 = Step2.active;
      actuators.V10 = false;
      actuators.V11 = Step3.active;
      actuators.V12 = Step5.active;
      actuators.V15 = Step8.active;
      actuators.V18 = Step10.active;
      actuators.V19 = false;
      actuators.V20 = Step13.active;
      actuators.V21 = false;
      actuators.V22 = Step10.active;
      actuators.V23 = Step10.active;
      actuators.V25 = Step13.active;
      actuators.V24 = Step13.active;

      connect(InitialStep1.outPort[1], Transition1.inPort) annotation (Line(
            points={{-159.5,100},{-144,100}}, color={0,0,0}));
      connect(Transition1.outPort, Step1.inPort[1]) annotation (Line(points={{
              -138.5,100},{-121,100}}, color={0,0,0}));
      connect(Step1.outPort[1], Transition2.inPort) 
        annotation (Line(points={{-99.5,100},{-84,100}}, color={0,0,0}));
      connect(Transition2.outPort, Step2.inPort[1]) 
        annotation (Line(points={{-78.5,100},{-61,100}}, color={0,0,0}));
      connect(Step2.outPort[1], Transition3.inPort) 
        annotation (Line(points={{-39.5,100},{-24,100}}, color={0,0,0}));
      connect(Transition3.outPort, Step3.inPort[1]) 
        annotation (Line(points={{-18.5,100},{-1,100}}, color={0,0,0}));
      connect(Step3.outPort[1], Transition4.inPort) 
        annotation (Line(points={{20.5,100},{36,100}}, color={0,0,0}));
      connect(Transition4.outPort, Step4.inPort[1]) 
        annotation (Line(points={{41.5,100},{59,100}}, color={0,0,0}));
      connect(Step4.outPort[1], Transition5.inPort) 
        annotation (Line(points={{80.5,100},{96,100}}, color={0,0,0}));
      connect(Transition5.outPort, Step5.inPort[1]) 
        annotation (Line(points={{101.5,100},{119,100}}, color={0,0,0}));
      connect(Step5.outPort[1], Transition6.inPort) 
        annotation (Line(points={{140.5,100},{156,100}}, color={0,0,0}));
      connect(Transition6.outPort, Step6.inPort[1]) annotation (Line(points={{
              161.5,100},{184,100},{184,70},{-160,70},{-160,40},{-121,40}},
            color={0,0,0}));
      connect(Step6.outPort[1], Transition7.inPort) 
        annotation (Line(points={{-99.5,40},{-84,40}}, color={0,0,0}));
      connect(Step12.inPort[1], Parallel1.split[1]) annotation (Line(points={{
              -63,-30},{-134.375,-30},{-134.375,-25}}, color={0,0,0}));
      connect(Step12.outPort[1], Transition12.inPort) 
        annotation (Line(points={{-41.5,-30},{-26,-30}}, color={0,0,0}));
      connect(Transition12.outPort, Step13.inPort[1]) 
        annotation (Line(points={{-20.5,-30},{-3,-30}}, color={0,0,0}));
      connect(Step13.outPort[1], Transition13.inPort) 
        annotation (Line(points={{18.5,-30},{34,-30}}, color={0,0,0}));
      connect(Transition13.outPort, Step14.inPort[1]) 
        annotation (Line(points={{39.5,-30},{57,-30}}, color={0,0,0}));
      connect(Step14.outPort[1], Parallel1.join[1]) annotation (Line(points={{
              78.5,-30},{152.375,-30},{152.375,-25}}, color={0,0,0}));
      connect(Step7.inPort[1], Parallel1.split[2]) annotation (Line(points={{
              -123,-70},{-138,-70},{-138,-75},{-134.375,-75}}, color={0,0,0}));
      connect(Step7.outPort[1], Transition8.inPort) annotation (Line(points={{
              -101.5,-70},{-86,-70}}, color={0,0,0}));
      connect(Transition8.outPort, Step8.inPort[1]) 
        annotation (Line(points={{-80.5,-70},{-63,-70}}, color={0,0,0}));
      connect(Step8.outPort[1], Transition9.inPort) 
        annotation (Line(points={{-41.5,-70},{-26,-70}}, color={0,0,0}));
      connect(Transition9.outPort, Step9.inPort[1]) 
        annotation (Line(points={{-20.5,-70},{-3,-70}}, color={0,0,0}));
      connect(Step9.outPort[1], Transition10.inPort) 
        annotation (Line(points={{18.5,-70},{34,-70}}, color={0,0,0}));
      connect(Transition10.outPort, Step10.inPort[1]) 
        annotation (Line(points={{39.5,-70},{57,-70}}, color={0,0,0}));
      connect(Step10.outPort[1], Transition11.inPort) 
        annotation (Line(points={{78.5,-70},{94,-70}}, color={0,0,0}));
      connect(Transition11.outPort, Step11.inPort[1]) 
        annotation (Line(points={{99.5,-70},{117,-70}}, color={0,0,0}));
      connect(Step11.outPort[1], Parallel1.join[2]) annotation (Line(points={{
              138.5,-70},{154,-70},{154,-75},{152.375,-75}}, color={0,0,0}));
      connect(Transition7.outPort, Parallel1.inPort) annotation (Line(points={{
              -78.5,40},{-40,40},{-40,20},{-190,20},{-190,-50},{-181.55,-50}},
            color={0,0,0}));
      connect(TransitionWithSignal1.inPort, Parallel1.outPort) annotation (Line(
            points={{2,-150},{208,-150},{208,-50},{197.7,-50}}, color={0,0,0}));
      connect(TransitionWithSignal1.outPort, InitialStep1.inPort[1]) annotation (Line(
            points={{-3.5,-150},{-194,-150},{-194,100},{-181,100}}, color={0,0,
              0}));
      connect(BooleanExpression1.y, TransitionWithSignal1.condition) annotation (Line(
            points={{-13.7,-132},{-2,-132},{-2,-138}}, color={255,0,255}));
    end Controller;

    package ControllerUtilities
      extends Modelica.Icons.Library;
      class Adapter_Inference
        Port_IdleTanks idleTanks;
      end Adapter_Inference;

      class Adapter_Superposition
        Port_Actuators actuators;
      end Adapter_Superposition;

      class Block_Recipe_TBD
        parameter Real startTime;
        parameter Real w_dilution=0.003;
        parameter Real w_concentrat=0.005;
        parameter Real T3_batch_level=0.1273;
        parameter Real T5_batch_level=0.211;
        Boolean trig;
        Boolean S0(start=true);
        Boolean S1;
        Boolean S2;
        Boolean S3;
        Boolean S4;
        Boolean S5;
        Boolean S6;
        Boolean S7;
        Boolean S8;
        Boolean S9;
        Boolean S10;
        Boolean S11;
        Boolean S12;
        Boolean S13;
        Boolean S14;
        Boolean tr0;
        Boolean tr1;
        Boolean tr2;
        Boolean tr3;
        Boolean tr4;
        Boolean tr5;
        Boolean tr6;
        Boolean tr7;
        Boolean tr8;
        Boolean tr9;
        Boolean tr10;
        Boolean tr11;
        Boolean tr12;
        Boolean tr13;
        Port_Actuators act annotation (Placement(transformation(extent={{-110,
                  -10},{-90,10}}, rotation=0)));
      end Block_Recipe_TBD;

      class BlockMain
        Boolean trig;

        Port_Actuators actuators annotation (Placement(transformation(extent={{
                  90,-10},{110,10}}, rotation=0)));
        Block_Recipe_TBD Recipe1 annotation (Placement(transformation(extent={{
                  -50,10},{-10,50}}, rotation=0)));
        Block_Recipe_TBD Recipe2 annotation (Placement(transformation(extent={{
                  10,10},{50,50}}, rotation=0)));
        Adapter_Inference Inference annotation (Placement(transformation(extent=
                 {{-50,-50},{-10,-10}}, rotation=0)));
        Adapter_Superposition Superposition annotation (Placement(
              transformation(extent={{10,-50},{50,-10}}, rotation=0)));
      end BlockMain;

      class Buffer_Recipe_TBD
        Port_Actuators act;
        Boolean S0;
        Boolean S1;
        Boolean S2;
        Boolean S3;
        Boolean S4;
        Boolean S5;
        Boolean S6;
        Boolean S7;
        Boolean S8;
        Boolean S9;
        Boolean S10;
        Boolean S11;
        Boolean S12;
        Boolean S13;
        Boolean S14;
      end Buffer_Recipe_TBD;

      class BufferMain
        Buffer_Recipe_TBD Recipe1;
        Buffer_Recipe_TBD Recipe2;
      end BufferMain;

      connector Port_Actuators
        output Boolean P1;
        output Boolean P2;
        output Boolean T5_Heater;
        output Boolean T7_Cooling;
        output Boolean T6_Cooling;
        output Boolean V1;
        output Boolean V2;
        output Boolean V3;
        output Boolean V4;
        output Boolean V5;
        output Boolean V6;
        output Boolean V8;
        output Boolean V9;
        output Boolean V10;
        output Boolean V11;
        output Boolean V12;
        output Boolean V15;
        output Boolean V18;
        output Boolean V19;
        output Boolean V20;
        output Boolean V21;
        output Boolean V22;
        output Boolean V23;
        output Boolean V24;
        output Boolean V25;

        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={Polygon(
                points={{-100,100},{100,0},{-100,-100},{-100,100}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}),
                                               Diagram(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics={Polygon(
                points={{0,50},{100,0},{0,-50},{0,50}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}));
      end Port_Actuators;

      connector Port_IdleTanks
        Boolean T5_idle;
        Boolean T7_idle;
      end Port_IdleTanks;

      connector Port_Sensors
        input Real LIS_301;
        input Real QI_302;
        input Real LIS_501;
        input Real QIS_502;
        input Real TI_503;
        input Real LIS_601;
        input Real TIS_602;
        input Real LIS_701;
        input Real TIS_702;
        annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={Polygon(
                points={{-100,100},{-100,-100},{100,0},{-100,100}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid)}),
                                               Diagram(coordinateSystem(
                preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics={Polygon(
                points={{0,50},{0,-50},{100,0},{0,50}},
                lineColor={0,0,0},
                lineThickness=0.5,
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid)}));
      end Port_Sensors;
    end ControllerUtilities;

    type Init = enumeration(
        GuessValues
          "GuessValues -- Guess values (not fixed) for p, T or h, X, C",
        InitialValues "InitialValues -- Initial values for p, T or h, X, C",
        SteadyStateMomentum "SteadyStateMomentum: Steady state momentum",
        SteadyStateHydraulic
          "SteadyStateHydraulic -- Hydraulic steady state (der(p)=0), guess value for p, initial values for T or h, X, C", 

        SteadyState
          "SteadyState -- Steady state (guess values for p, T or h, X, C)")
      "Enumeration to define initialization options" 
    annotation (Documentation(info="<html>
<p>
Integer type that can have the following values
(to be selected via choices menu):
</p>
 
<table border=1 cellspacing=0 cellpadding=2>
<tr><th><b>Types.Init.</b></th><th><b>Meaning</b></th></tr>
<tr><td>GuessValues</td>
    <td>GuessValues -- Guess values (not fixed) for p, T or h, X, C</td></tr>
 
<tr><td>InitialValues</td>
    <td>Initial values for p, T or h, X, C</td></tr>
 
<tr><td>SteadyStateMomentum</td>
    <td>Steady state momentum</td></tr>
 
<tr><td>SteadyStateHydraulic</td>
    <td>Hydraulic steady state (der(p)=0), guess value for p, 
        initial values for T or h, X, C</td></tr>
 
<tr><td>SteadyState</td>
    <td>Steady state (guess values for p, T or h, X, C)</td></tr>
</table>
</html>"));

  end BaseClasses;

end AST_BatchPlant;
