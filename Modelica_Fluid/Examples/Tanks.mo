within Modelica_Fluid.Examples;
package Tanks "Library demonstrating the usage of the tank model"
  extends Modelica.Icons.Library;
  model OneTank
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    Modelica_Fluid.Volumes.Tank tank(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      area=1,
      levelMax=1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)},
      level_start=0.9,
      V0=0.1,
      nTopPorts=1,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{-40,28},{0,68}}, rotation=0)));

    Sources.PrescribedMassFlowRate_TX flowSource(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.OutOfPort,
      m_flow=20,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{-52,70},{-32,90}}, rotation
            =0)));
    annotation (Diagram(graphics),
      experiment(StopTime=100),
      experimentSetupOutput,
      Commands(file="../Scripts/Examples/OneTank/plot level and port.p.mos"
          "plot level and port.p", file=
            "../Scripts/Examples/OneTank/plot level, port.p and port.m_flow.mos"
          "plot level, port.p and port.m_flow"));
    inner Ambient ambient annotation (Placement(transformation(extent={{-10,72},
              {10,92}}, rotation=0)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient_fixed(
                                           redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.InToPort,
      p=ambient.default_p_ambient,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{-54,-20},{-34,0}}, rotation
            =0)));
    PressureLosses.WallFrictionAndGravity pipe(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={-20,10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(flowSource.port, tank.topPorts[1])  annotation (Line(points={{-32,
            80},{-20,80},{-20,69}}, color={0,127,255}));
    connect(ambient_fixed.port, pipe.port_a) annotation (Line(points={{-34,-10},
            {-20,-10},{-20,0}}, color={0,127,255}));
    connect(pipe.port_b, tank.ports[1]) annotation (Line(points={{-20,20},{-20,
            27}}, color={0,127,255}));
  end OneTank;

  model TwoTanks
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;
    parameter Boolean stiffCharacteristicForEmptyPort=true;

    annotation (Diagram(graphics),
      experiment(StopTime=70),
      experimentSetupOutput,
      Commands(file="../Scripts/Examples/TwoTanks/plot level and port.p.mos"
          "plot level and port.p"));
    inner Ambient ambient annotation (Placement(transformation(extent={{40,62},
              {60,82}}, rotation=0)));
    Modelica_Fluid.Volumes.Tank tank1(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort,
      area=1,
      levelMax=4,
      level_start=3,
      T_start=Modelica.SIunits.Conversions.from_degC(50),
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)}) 
      annotation (Placement(transformation(extent={{-80,0},{-40,40}}, rotation=
              0)));
    Modelica_Fluid.Volumes.Tank tank2(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort,
      area=1,
      levelMax=4,
      level_start=1,
      T_start=Modelica.SIunits.Conversions.from_degC(100),
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)}) 
      annotation (Placement(transformation(extent={{0,0},{40,40}}, rotation=0)));
    PressureLosses.WallFrictionAndGravity pipe(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1)  annotation (Placement(transformation(extent={{-30,-30},{
              -10,-10}}, rotation=0)));
  equation
    connect(tank1.ports[1], pipe.port_a) annotation (Line(points={{-60,-1},{-60,
            -20},{-30,-20}}, color={0,127,255}));
    connect(pipe.port_b, tank2.ports[1]) annotation (Line(points={{-10,-20},{20,
            -20},{20,-1}}, color={0,127,255}));
  end TwoTanks;

  model TankWithEmptyingPipe1
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    Sources.PrescribedMassFlowRate_TX flowSource(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.OutOfPort,
      m_flow=50,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{-20,40},{0,60}}, rotation=0)));
    annotation (Diagram(graphics),
      experiment(StopTime=35),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/TankWithEmptyingPipe1/plot level and port.p.mos"
          "plot level and port.p", file=
            "../Scripts/Examples/TankWithEmptyingPipe1/plot level, port.p and port.m_flow.mos"
          "plot level, port.p and port.m_flow"));
    inner Ambient ambient annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient_fixed(
                                           redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.InToPort,
      p=ambient.default_p_ambient,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{-60,-100},{-40,-80}},
            rotation=0)));
    ControlValves.ValveDiscrete valveDiscrete(redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater, Kv=100) 
      annotation (Placement(transformation(
          origin={-20,-50},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica.Blocks.Sources.BooleanConstant open(k=false) 
      annotation (Placement(transformation(extent={{-60,-60},{-40,-40}},
            rotation=0)));
    Modelica_Fluid.Volumes.Tank tank1(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      area=1,
      V0=0.1,
      levelMax=2,
      level_start=0.1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.1,
          portLevel=1)},
      stiffCharacteristicForEmptyPort=true,
      nTopPorts=1) 
      annotation (Placement(transformation(extent={{-40,-20},{0,20}}, rotation=
              0)));
    PressureLosses.WallFrictionAndGravity pipe(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={40,10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(ambient_fixed.port, valveDiscrete.port_a) annotation (Line(points={
            {-40,-90},{-20,-90},{-20,-60}}, color={0,127,255}));
    connect(open.y, valveDiscrete.open) annotation (Line(points={{-39,-50},{-28,
            -50}}, color={255,0,255}));
    connect(flowSource.port, pipe.port_b) annotation (Line(points={{0,50},{40,
            50},{40,20}}, color={0,127,255}));
    connect(valveDiscrete.port_b, tank1.ports[1]) annotation (Line(points={{-20,
            -40},{-20,-21}}, color={0,127,255}));
    connect(pipe.port_a, tank1.ports[2]) annotation (Line(points={{40,0},{40,
            -28},{-18,-28},{-18,-20},{-20,-20},{-20,-21}}, color={0,127,255}));
  end TankWithEmptyingPipe1;

  model TankWithEmptyingPipe2
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    annotation (Diagram(graphics),
      experiment(StopTime=35),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/TankWithEmptyingPipe2/plot level and port.p.mos"
          "plot level and port.p"));
    inner Ambient ambient annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient_fixed(
                                           redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.InToPort,
      p=ambient.default_p_ambient,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{-60,-100},{-40,-80}},
            rotation=0)));
    Modelica_Fluid.Volumes.Tank tank1(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      area=1,
      V0=0.1,
      levelMax=2,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.1,
          portLevel=1)},
      level_start=2,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{-40,-20},{0,20}}, rotation=
              0)));
    PressureLosses.WallFrictionAndGravity pipe1(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={-20,-60},
          extent={{-10,-10},{10,10}},
          rotation=90)));

    PressureLosses.WallFrictionAndGravity pipe2(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={30,-60},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient_fixed1(
                                           redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.InToPort,
      p=ambient.default_p_ambient,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{0,-100},{20,-80}}, rotation
            =0)));
  equation
    connect(tank1.ports[1], pipe1.port_b) annotation (Line(points={{-20,-21},{
            -20,-50}}, color={0,127,255}));
    connect(ambient_fixed.port, pipe1.port_a) annotation (Line(points={{-40,-90},
            {-20,-90},{-20,-70}}, color={0,127,255}));
    connect(tank1.ports[2], pipe2.port_b) annotation (Line(points={{-20,-21},{
            -18,-21},{-18,-40},{30,-40},{30,-50}}, color={0,127,255}));
    connect(ambient_fixed1.port, pipe2.port_a) annotation (Line(points={{20,-90},
            {30,-90},{30,-70}}, color={0,127,255}));
  end TankWithEmptyingPipe2;

  model TanksWithEmptyingPipe1
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    annotation (Diagram(graphics),
      experiment(StopTime=35),
      experimentSetupOutput,
      Commands(
        file=
            "../Scripts/Examples/TanksWithEmptyingPipe1/plot level, port.p and port.m_flow.mos"
          "plot level, port.p and port.m_flow"));
    inner Ambient ambient annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient_fixed1(
                                            redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater, flowDirection=
          Modelica_Fluid.Types.SourceFlowDirection.InToPort,
      p=ambient.default_p_ambient,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{-100,-80},{-80,-60}},
            rotation=0)));
    Modelica_Fluid.Volumes.Tank tank1(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      area=1,
      V0=0.1,
      levelMax=2,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.1,
          portLevel=1)},
      level_start=2,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{-80,0},{-40,40}}, rotation=
              0)));
    PressureLosses.WallFrictionAndGravity pipe1(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={-60,-40},
          extent={{-10,-10},{10,10}},
          rotation=90)));

    PressureLosses.WallFrictionAndGravity pipe2(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={40,-40},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient_fixed2(
                                           redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.InToPort,
      p=ambient.default_p_ambient,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{0,-80},{20,-60}}, rotation=
              0)));
    Modelica_Fluid.Volumes.Tank tank2(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      area=1,
      V0=0.1,
      levelMax=2,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.1,
          portLevel=0.5)},
      level_start=0.1,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{20,0},{60,40}}, rotation=0)));
    PressureLosses.WallFrictionAndGravity pipe3(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed,
      length=1,
      diameter=0.1,
      height_ab=-0.5) 
                   annotation (Placement(transformation(extent={{-20,10},{0,30}}, 
            rotation=0)));
  equation
    connect(tank1.ports[1], pipe1.port_b) annotation (Line(points={{-60,-1},{
            -60,-30}}, color={0,127,255}));
    connect(ambient_fixed1.port, pipe1.port_a) 
                                              annotation (Line(points={{-80,-70},
            {-60,-70},{-60,-50}}, color={0,127,255}));
    connect(ambient_fixed2.port, pipe2.port_a) annotation (Line(points={{20,-70},
            {40,-70},{40,-50}}, color={0,127,255}));
    connect(tank2.ports[1], pipe2.port_b) 
      annotation (Line(points={{40,-1},{40,-30}}, color={0,127,255}));
    connect(pipe3.port_a, tank1.ports[2]) annotation (Line(points={{-20,20},{
            -30,20},{-30,-10},{-58,-10},{-58,0},{-60,0},{-60,-1}}, color={0,127,
            255}));
    connect(pipe3.port_b, tank2.ports[2]) annotation (Line(points={{0,20},{10,
            20},{10,-8},{38,-8},{38,0},{40,0},{40,-1}}, color={0,127,255}));
  end TanksWithEmptyingPipe1;

  model TanksWithEmptyingPipe2
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    parameter Boolean stiffCharacteristicForEmptyPort=true;
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;
    replaceable package Medium = 
       Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater                    constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    annotation (Diagram(graphics),
      experiment(StopTime=300),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/TanksWithEmptyingPipe2/plot level and port.m_flow.mos"
          "plot level and port.m_flow"));
    inner Ambient ambient annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.FixedBoundary_pTX ambient_fixed(
                                           redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      flowDirection=Modelica_Fluid.Types.SourceFlowDirection.InToPort,
      p=ambient.default_p_ambient,
      T=ambient.default_T_ambient) 
      annotation (Placement(transformation(extent={{-16,-102},{-36,-82}},
            rotation=0)));
    ControlValves.ValveDiscrete valveDiscrete(redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater, Kv=100) 
      annotation (Placement(transformation(
          origin={-60,-78},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica.Blocks.Sources.BooleanConstant open(k=false) 
      annotation (Placement(transformation(extent={{-98,-88},{-78,-68}},
            rotation=0)));
    Modelica_Fluid.Volumes.Tank tank3(
      redeclare package Medium = Medium,
      area=1,
      V0=0.1,
      levelMax=20,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.05,
          portLevel=6.5)},
      level_start=6,
      nTopPorts=1,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{-80,-50},{-40,-10}},
            rotation=0)));
    Modelica_Fluid.Volumes.Tank tank1(
      redeclare package Medium = Medium,
      area=1,
      V0=0.1,
      levelMax=10,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)},
      level_start=9,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{50,50},{90,90}}, rotation=0)));
    Modelica_Fluid.Volumes.Tank tank2(
      redeclare package Medium = Medium,
      area=1,
      V0=0.1,
      levelMax=10,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Volumes.BaseClasses.TankPortData(diameter=0.05,
          portLevel=2),Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.1, portLevel=3)},
      level_start=1,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{-20,10},{20,50}}, rotation=
              0)));
    PressureLosses.StaticHead pipe1(redeclare package Medium = Medium,
        height_ab=2) annotation (Placement(transformation(
          origin={70,30},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    PressureLosses.StaticHead pipe2(redeclare package Medium = Medium,
        height_ab=2) annotation (Placement(transformation(
          origin={0,-24},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    PressureLosses.StaticHead pipe3(redeclare package Medium = Medium,
        height_ab=2) annotation (Placement(transformation(
          origin={-60,10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(ambient_fixed.port, valveDiscrete.port_a) annotation (Line(points={
            {-36,-92},{-60,-92},{-60,-88}}, color={0,127,255}));
    connect(open.y, valveDiscrete.open) annotation (Line(points={{-77,-78},{-68,
            -78}}, color={255,0,255}));
    connect(valveDiscrete.port_b,tank3. ports[1]) annotation (Line(points={{-60,
            -68},{-60,-51}}, color={0,127,255}));
    connect(pipe1.port_b, tank1.ports[1]) annotation (Line(points={{70,40},{70,
            49}}, color={0,127,255}));
    connect(pipe2.port_b, tank2.ports[2]) annotation (Line(points={{
            6.12323e-016,-14},{6.12323e-016,-10},{0,-10},{0,9}}, color={0,127,
            255}));
    connect(pipe2.port_a, tank3.ports[2]) annotation (Line(points={{
            -6.12323e-016,-34},{-6.12323e-016,-48},{0,-60},{-58,-60},{-58,-51},
            {-60,-51}}, color={0,127,255}));
    connect(pipe3.port_a, tank3.topPorts[1]) 
                                            annotation (Line(points={{-60,0},{
            -60,-9}}, color={0,127,255}));
    connect(pipe3.port_b, tank2.ports[1]) annotation (Line(points={{-60,20},{
            -60,26},{-30,26},{-30,0},{-2,0},{-2,9},{0,9}}, color={0,127,255}));
    connect(pipe1.port_a, tank2.ports[3]) annotation (Line(points={{70,20},{70,
            0},{2,0},{2,9},{0,9}}, color={0,127,255}));
  end TanksWithEmptyingPipe2;

  model ThreeOpenTanks "Demonstrating the usage of OpenTank"
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
     // replaceable package Medium = Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater extends
    // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends
    // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
     replaceable package Medium = 
        Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater                           constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    Modelica_Fluid.Volumes.OpenTank tank1(
      area=1,
      redeclare package Medium = Medium,
      p_static_at_port=false,
      height=12,
      level_start=8,
      zeta_in={1.05},
      n_ports=1,
      pipe_diameters={0.1}) 
                     annotation (Placement(transformation(extent={{-80,20},{-40,
              60}}, rotation=0)));
    Modelica_Fluid.Volumes.OpenTank tank2(
      area=1,
      redeclare package Medium = Medium,
      p_static_at_port=false,
      height=12,
      level_start=3,
      zeta_in={1.05},
      n_ports=1,
      pipe_diameters={0.1}) 
                     annotation (Placement(transformation(extent={{-20,20},{20,
              60}}, rotation=0)));
    annotation (Diagram(graphics),
      experiment(StopTime=100),
      experimentSetupOutput,
      Documentation(info="<html>
  
</html>"));

    inner Modelica_Fluid.Ambient ambient 
                                     annotation (Placement(transformation(
            extent={{76,-96},{96,-76}}, rotation=0)));
    Modelica_Fluid.Volumes.OpenTank tank3(
      area=1,
      redeclare package Medium = Medium,
      p_static_at_port=false,
      height=12,
      level_start=3,
      zeta_in={1.05},
      n_ports=1,
      pipe_diameters={0.1}) 
                     annotation (Placement(transformation(extent={{40,20},{80,
              60}}, rotation=0)));
    Modelica_Fluid.PressureLosses.StaticHead pipe1(           redeclare package
        Medium =                                                                       Medium,
      flowDirection=Modelica_Fluid.Types.FlowDirection.Bidirectional,
      height_ab=2) annotation (Placement(transformation(
          origin={-60,-10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.PressureLosses.StaticHead pipe2(           redeclare package
        Medium =                                                                       Medium,
      flowDirection=Modelica_Fluid.Types.FlowDirection.Bidirectional,
      height_ab=2) annotation (Placement(transformation(
          origin={0,-10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.PressureLosses.StaticHead pipe3(           redeclare package
        Medium =                                                                       Medium,
      flowDirection=Modelica_Fluid.Types.FlowDirection.Bidirectional,
      height_ab=-1) annotation (Placement(transformation(
          origin={60,-10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(pipe1.port_a, pipe2.port_a) annotation (Line(points={{-60,-20},{-60,
            -40},{-6.12323e-016,-40},{-6.12323e-016,-20}}, color={0,127,255}));
    connect(pipe2.port_a, pipe3.port_a) annotation (Line(points={{-6.12323e-016,
            -20},{0,-20},{0,-40},{60,-40},{60,-20}}, color={0,127,255}));
    connect(pipe2.port_b, tank2.ports[1]) annotation (Line(points={{
            6.12323e-016,0},{6.12323e-016,10},{0,10},{0,19}}, color={0,127,255}));
    connect(pipe3.port_b, tank3.ports[1]) 
      annotation (Line(points={{60,0},{60,19}}, color={0,127,255}));
    connect(pipe1.port_b, tank1.ports[1]) annotation (Line(points={{-60,0},{-60,
            19}}, color={0,127,255}));
  end ThreeOpenTanks;

  model TestEmptyOpenTank "Test whether an empty tank is properly handeled"
    extends Modelica.Icons.Example;
    Modelica_Fluid.Volumes.Tank tank1(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      area=1,
      level_start=1,
      levelMax=1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)}) 
                              annotation (Placement(transformation(extent={{-20,
              20},{20,60}}, rotation=0)));

    Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe(
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      redeclare package WallFriction = 
          Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.LaminarAndQuadraticTurbulent,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={0,-10},
          extent={{10,-10},{-10,10}},
          rotation=270)));

    annotation (
      Diagram(graphics),
      experiment(StopTime=25),
      experimentSetupOutput);
    Modelica_Fluid.Volumes.Tank tank2(
      area=1,
      level_start=0,
      redeclare package Medium = 
          Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater,
      levelMax=1,
      portsData={Modelica_Fluid.Volumes.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)},
      initType=Modelica_Fluid.Types.Init.NoInit) 
      annotation (Placement(transformation(extent={{-20,-80},{20,-40}},
            rotation=0)));
    inner Modelica_Fluid.Ambient ambient 
                                     annotation (Placement(transformation(
            extent={{56,58},{76,78}}, rotation=0)));
  equation
    connect(pipe.port_b, tank1.ports[1]) annotation (Line(points={{6.12323e-016,
            0},{0,0},{0,19}}, color={0,127,255}));
    connect(pipe.port_a, tank2.topPorts[1]) annotation (Line(points={{
            -6.12323e-016,-20},{0,-20},{0,-39}}, color={0,127,255}));
  end TestEmptyOpenTank;

end Tanks;
