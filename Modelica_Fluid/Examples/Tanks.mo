within Modelica_Fluid.Examples;
package Tanks "Library demonstrating the usage of the tank model"
  extends Modelica.Icons.Library;
  model OneTank
    "Tank with one time-varying top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    Modelica_Fluid.Vessels.TankWithTopPorts tank(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      height=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.1, height=0)},
      V0=0.1,
      nTopPorts=1,
      nPorts=1,
      level_start=0) 
      annotation (Placement(transformation(extent={{0,0},{40,40}},   rotation=0)));

    Sources.MassFlowSource_T flowSource(nPorts=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      m_flow=20,
      T=system.T_ambient,
      use_m_flow_in=true) 
      annotation (Placement(transformation(extent={{-12,42},{8,62}},   rotation=
             0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=100),
      experimentSetupOutput,
      Commands(file="../Scripts/Examples/OneTank/plot level and port.p.mos"
          "plot level and port.p", file=
            "../Scripts/Examples/OneTank/plot level, port.p and port.m_flow.mos"
          "plot level, port.p and port.m_flow"));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{70,72},
              {90,92}}, rotation=0)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed(nPorts=1,
                                           redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-14,-50},{6,-30}}, rotation=
             0)));
    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=-1) 
                   annotation (Placement(transformation(
          origin={20,-18},
          extent={{10,-10},{-10,10}},
          rotation=90)));
    Modelica.Blocks.Sources.TimeTable timeTable(table=[0,0; 10,0; 10,40; 20,40;
          20,10; 50,10; 50,0; 60,0; 60,20; 70,20; 80,55; 80,0; 100,0]) 
      annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
  equation
    connect(flowSource.ports[1], tank.topPorts[1])  annotation (Line(points={{8,52},{
            20,52},{20,41}},        color={0,127,255}));
    connect(tank.ports[1], pipe.port_a) annotation (Line(
        points={{20,-1},{20,-8}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_b, ambient_fixed.ports[1]) annotation (Line(
        points={{20,-28},{20,-40},{6,-40}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(timeTable.y, flowSource.m_flow_in) annotation (Line(
        points={{-39,70},{-24,70},{-24,60},{-12,60}},
        color={0,0,127},
        smooth=Smooth.None));
  end OneTank;

  model TwoTanks
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;
    parameter Boolean stiffCharacteristicForEmptyPort=true;

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=70),
      experimentSetupOutput,
      Commands(file="../Scripts/Examples/TwoTanks/plot level and port.p.mos"
          "plot level and port.p"));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{40,62},
              {60,82}}, rotation=0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort,
      crossArea=1,
      height=4,
      level_start=3,
      T_start=Modelica.SIunits.Conversions.from_degC(50),
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.1, height=0)}) 
      annotation (Placement(transformation(extent={{-80,0},{-40,40}}, rotation=
              0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank2(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort,
      crossArea=1,
      height=4,
      level_start=1,
      T_start=Modelica.SIunits.Conversions.from_degC(100),
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.1, height=0)}) 
      annotation (Placement(transformation(extent={{0,0},{40,40}}, rotation=0)));
    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
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

    Sources.MassFlowSource_T flowSource(
      nPorts=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      m_flow=50,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-20,40},{0,60}}, rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=35),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/TankWithEmptyingPipe1/plot level and port.p.mos"
          "plot level and port.p", file=
            "../Scripts/Examples/TankWithEmptyingPipe1/plot level, port.p and port.m_flow.mos"
          "plot level, port.p and port.m_flow"));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed(nPorts=1,
                                           redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-60,-100},{-40,-80}},
            rotation=0)));
    Modelica_Fluid.Valves.ValveDiscrete valveDiscrete(
                                              redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      dp_nominal(displayUnit="Pa") = 1,
      m_flow_nominal=100) 
      annotation (Placement(transformation(
          origin={-20,-50},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica.Blocks.Sources.BooleanConstant open(k=false) 
      annotation (Placement(transformation(extent={{-60,-60},{-40,-40}},
            rotation=0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      level_start=0.1,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.05, height=0),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
                                                          diameter=0.1,
          height=1)},
      stiffCharacteristicForEmptyPort=true,
      nTopPorts=1) 
      annotation (Placement(transformation(extent={{-40,-20},{0,20}}, rotation=
              0)));
    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={40,10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(ambient_fixed.ports[1], valveDiscrete.port_a) annotation (Line(points={
            {-40,-90},{-20,-90},{-20,-60}}, color={0,127,255}));
    connect(open.y, valveDiscrete.open) annotation (Line(points={{-39,-50},{-28,
            -50}}, color={255,0,255}));
    connect(flowSource.ports[1], pipe.port_b) annotation (Line(points={{0,50},{40,
            50},{40,20}}, color={0,127,255}));
    connect(valveDiscrete.port_b, tank1.ports[1]) annotation (Line(points={{-20,-40},
            {-20,-30},{-20,-21},{-22,-21}},
                             color={0,127,255}));
    connect(pipe.port_a, tank1.ports[2]) annotation (Line(points={{40,0},{40,
            -28},{-18,-28},{-18,-20},{-18,-20},{-18,-21}}, color={0,127,255}));
  end TankWithEmptyingPipe1;

  model TankWithEmptyingPipe2
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=35),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/TankWithEmptyingPipe2/plot level and port.p.mos"
          "plot level and port.p"));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed(nPorts=1,
                                           redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-60,-100},{-40,-80}},
            rotation=0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.05, height=0),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
                                                          diameter=0.1,
          height=1)},
      level_start=2,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{-40,-20},{0,20}}, rotation=
              0)));
    Modelica_Fluid.Pipes.StaticPipe pipe1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={-20,-60},
          extent={{-10,-10},{10,10}},
          rotation=90)));

    Modelica_Fluid.Pipes.StaticPipe pipe2(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={30,-60},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed1(nPorts=1,
                                           redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{0,-100},{20,-80}}, rotation=
             0)));
  equation
    connect(tank1.ports[1], pipe1.port_b) annotation (Line(points={{-22,-21},{
            -22,-35},{-20,-35},{-20,-50}},
                       color={0,127,255}));
    connect(ambient_fixed.ports[1], pipe1.port_a) annotation (Line(points={{-40,-90},
            {-20,-90},{-20,-70}}, color={0,127,255}));
    connect(tank1.ports[2], pipe2.port_b) annotation (Line(points={{-18,-21},{
            -18,-21},{-18,-40},{30,-40},{30,-50}}, color={0,127,255}));
    connect(ambient_fixed1.ports[1], pipe2.port_a) annotation (Line(points={{20,-90},
            {30,-90},{30,-70}}, color={0,127,255}));
  end TankWithEmptyingPipe2;

  model TanksWithEmptyingPipe1
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=35),
      experimentSetupOutput,
      Commands(
        file=
            "../Scripts/Examples/TanksWithEmptyingPipe1/plot level, port.p and port.m_flow.mos"
          "plot level, port.p and port.m_flow"));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed1(nPorts=1,
                                            redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-100,-80},{-80,-60}},
            rotation=0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.05, height=0),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
                                                          diameter=0.1,
          height=1)},
      level_start=2,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{-80,0},{-40,40}}, rotation=
              0)));
    Modelica_Fluid.Pipes.StaticPipe pipe1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={-60,-40},
          extent={{-10,-10},{10,10}},
          rotation=90)));

    Modelica_Fluid.Pipes.StaticPipe pipe2(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={40,-40},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed2(nPorts=1,
                                           redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{0,-80},{20,-60}}, rotation=
              0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank2(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.05, height=0),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
                                                          diameter=0.1,
          height=0.5)},
      level_start=0.1,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{20,0},{60,40}}, rotation=0)));
    Modelica_Fluid.Pipes.StaticPipe pipe3(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=-0.5) 
                   annotation (Placement(transformation(extent={{-20,10},{0,30}},
            rotation=0)));
  equation
    connect(tank1.ports[1], pipe1.port_b) annotation (Line(points={{-62,-1},{
            -62,-15},{-60,-15},{-60,-30}},
                       color={0,127,255}));
    connect(ambient_fixed1.ports[1], pipe1.port_a) 
                                              annotation (Line(points={{-80,-70},
            {-60,-70},{-60,-50}}, color={0,127,255}));
    connect(ambient_fixed2.ports[1], pipe2.port_a) annotation (Line(points={{20,-70},
            {40,-70},{40,-50}}, color={0,127,255}));
    connect(tank2.ports[1], pipe2.port_b) 
      annotation (Line(points={{38,-1},{38,-15},{40,-15},{40,-30}},
                                                  color={0,127,255}));
    connect(pipe3.port_a, tank1.ports[2]) annotation (Line(points={{-20,20},{
            -30,20},{-30,-10},{-58,-10},{-58,0},{-58,0},{-58,-1}}, color={0,127,
            255}));
    connect(pipe3.port_b, tank2.ports[2]) annotation (Line(points={{0,20},{10,
            20},{10,-8},{38,-8},{38,0},{42,0},{42,-1}}, color={0,127,255}));
  end TanksWithEmptyingPipe1;

  model TanksWithEmptyingPipe2
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    parameter Boolean stiffCharacteristicForEmptyPort=true;
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;
    replaceable package Medium = 
       Modelica.Media.Water.ConstantPropertyLiquidWater                    constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=300),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/TanksWithEmptyingPipe2/plot level and port.m_flow.mos"
          "plot level and port.m_flow"));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-100,60},
              {-80,80}}, rotation=0)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed(nPorts=1,
                                           redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-16,-102},{-36,-82}},
            rotation=0)));
    Modelica_Fluid.Valves.ValveDiscrete valveDiscrete(
                                              redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      dp_nominal(displayUnit="Pa") = 1,
      m_flow_nominal=100) 
      annotation (Placement(transformation(
          origin={-60,-78},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica.Blocks.Sources.BooleanConstant open(k=false) 
      annotation (Placement(transformation(extent={{-98,-88},{-78,-68}},
            rotation=0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank3(
      redeclare package Medium = Medium,
      crossArea=1,
      V0=0.1,
      height=20,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.05, height=0),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
                                                          diameter=0.05,
          height=6.5)},
      level_start=6,
      nTopPorts=1,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{-80,-50},{-40,-10}},
            rotation=0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank1(
      redeclare package Medium = Medium,
      crossArea=1,
      V0=0.1,
      height=10,
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.1, height=0)},
      level_start=9,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{50,50},{90,90}}, rotation=0)));
    Modelica_Fluid.Vessels.TankWithTopPorts tank2(
      redeclare package Medium = Medium,
      crossArea=1,
      V0=0.1,
      height=10,
      nPorts=3,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.05, height=0),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
                                                          diameter=0.05,
          height=2),Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(
          diameter=0.1, height=3)},
      level_start=1,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{-20,10},{20,50}}, rotation=
              0)));
    Modelica_Fluid.Fittings.GenericStaticHead pipe1(
                                    redeclare package Medium = Medium,
        height_ab=2) annotation (Placement(transformation(
          origin={70,30},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Fittings.GenericStaticHead pipe2(
                                    redeclare package Medium = Medium,
        height_ab=2) annotation (Placement(transformation(
          origin={0,-22},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Fittings.GenericStaticHead pipe3(
                                    redeclare package Medium = Medium,
        height_ab=2) annotation (Placement(transformation(
          origin={-60,10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(ambient_fixed.ports[1], valveDiscrete.port_a) annotation (Line(points={
            {-36,-92},{-60,-92},{-60,-88}}, color={0,127,255}));
    connect(open.y, valveDiscrete.open) annotation (Line(points={{-77,-78},{-68,
            -78}}, color={255,0,255}));
    connect(valveDiscrete.port_b,tank3. ports[1]) annotation (Line(points={{-60,-68},
            {-60,-59},{-60,-51},{-62,-51}},
                             color={0,127,255}));
    connect(pipe1.port_b, tank1.ports[1]) annotation (Line(points={{70,40},{70,
            45},{70,49},{70,49}},
                  color={0,127,255}));
    connect(pipe2.port_a, tank3.ports[2]) annotation (Line(points={{
            -6.12323e-016,-32},{-6.12323e-016,-48},{0,-60},{-58,-60},{-58,-51},
            {-58,-51}}, color={0,127,255}));
    connect(pipe3.port_a, tank3.topPorts[1]) 
                                            annotation (Line(points={{-60,0},{
            -60,-5},{-60,-9},{-60,-9}},
                      color={0,127,255}));
    connect(pipe3.port_b, tank2.ports[1]) annotation (Line(points={{-60,20},{
            -60,26},{-30,26},{-30,0},{-2,0},{-2,9},{-2.66667,9}},
                                                           color={0,127,255}));
    connect(pipe1.port_a, tank2.ports[3]) annotation (Line(points={{70,20},{70,
            0},{2,0},{2,9},{2.66667,9}},
                                   color={0,127,255}));
    connect(pipe2.port_b, tank2.ports[2]) annotation (Line(
        points={{6.12323e-016,-12},{0,-12},{0,9},{2.22045e-016,9}},
        color={0,127,255},
        smooth=Smooth.None));
  end TanksWithEmptyingPipe2;

  model ThreeTanks "Demonstrating the usage of SimpleTank"
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
     // replaceable package Medium = Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater extends
    // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends
    // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
     replaceable package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater                           constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    Modelica_Fluid.Vessels.SimpleTank tank1(
      crossArea=1,
      redeclare package Medium = Medium,
      use_portsData=true,
      height=12,
      level_start=8,
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=
          0.1)})     annotation (Placement(transformation(extent={{-80,20},{-40,
              60}}, rotation=0)));
    Modelica_Fluid.Vessels.SimpleTank tank2(
      crossArea=1,
      redeclare package Medium = Medium,
      use_portsData=true,
      height=12,
      level_start=3,
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=
          0.1)})     annotation (Placement(transformation(extent={{-20,20},{20,
              60}}, rotation=0)));

    inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.FixedInitial) 
                                     annotation (Placement(transformation(
            extent={{70,-90},{90,-70}}, rotation=0)));
    Modelica_Fluid.Vessels.SimpleTank tank3(
      crossArea=1,
      redeclare package Medium = Medium,
      use_portsData=true,
      height=12,
      level_start=3,
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=
          0.1)})     annotation (Placement(transformation(extent={{40,10},{80,50}},
                    rotation=0)));
    Modelica_Fluid.Fittings.GenericStaticHead pipe1(          redeclare package
        Medium =                                                                       Medium,
      allowFlowReversal=true,
      height_ab=2) annotation (Placement(transformation(
          origin={-60,-10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Fittings.GenericStaticHead pipe2(          redeclare package
        Medium =                                                                       Medium,
      allowFlowReversal=true,
      height_ab=2) annotation (Placement(transformation(
          origin={0,-10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
    Modelica_Fluid.Fittings.GenericStaticHead pipe3(          redeclare package
        Medium =                                                                       Medium,
      allowFlowReversal=true,
      height_ab=-1) annotation (Placement(transformation(
          origin={60,-20},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(pipe1.port_a, pipe2.port_a) annotation (Line(points={{-60,-20},{-60,
            -40},{-6.12323e-016,-40},{-6.12323e-016,-20}}, color={0,127,255}));
    connect(pipe2.port_a, pipe3.port_a) annotation (Line(points={{-6.12323e-016,
            -20},{0,-20},{0,-40},{60,-40},{60,-30}}, color={0,127,255}));
    connect(pipe3.port_b, tank3.ports[1]) 
      annotation (Line(points={{60,-10},{60,-10},{60,10}},
                                                color={0,127,255}));
    connect(pipe1.port_b, tank1.ports[1]) annotation (Line(points={{-60,0},{-60,
            10},{-60,20}},
                  color={0,127,255}));
    connect(pipe2.port_b, tank2.ports[1]) annotation (Line(
        points={{6.12323e-016,0},{0,0},{0,20}},
        color={0,127,255},
        smooth=Smooth.None));

    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=200),
      experimentSetupOutput,
      Commands(file=
            "../Scripts/Examples/ThreeTanks/plot level and port.m_flow.mos"
          "plot level and port.m_flow"),
      Documentation(info="<html>
  
</html>"));

  end ThreeTanks;

  model TanksWithOverflow "Two tanks connected with pipes at different heights"
    extends Modelica.Icons.Example;
    import Modelica_Fluid;
    Modelica_Fluid.Vessels.SimpleTank upperTank(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      height=20,
      level_start=2,
      crossArea=0.2,
      nPorts=3,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1),
          Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1, height=
           10)}) 
      annotation (Placement(transformation(extent={{-40,20},{0,60}}, rotation=0)));
    Modelica_Fluid.Sources.MassFlowSource_T massFlowRate(nPorts=1,
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      m_flow=0.2,
      use_m_flow_in=true) 
      annotation (Placement(transformation(extent={{-60,-40},{-40,-20}}, rotation=
             0)));
    inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.FixedInitial) 
                                        annotation (Placement(transformation(
            extent={{-150,-112},{-130,-92}},  rotation=0)));
    Modelica_Fluid.Sensors.Pressure pressure(redeclare package Medium = 
          Modelica.Media.Water.StandardWater) annotation (Placement(
          transformation(extent={{40,16},{60,36}}, rotation=0)));
    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      diameter=0.02,
      height_ab=-20,
      length=200) annotation (Placement(transformation(
          origin={0,-30},
          extent={{10,-10},{-10,10}},
          rotation=90)));

    Modelica_Fluid.Vessels.SimpleTank lowerTank(
      height=20,
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      level_start=2,
      crossArea=1,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=
          0.1),Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=0.1,
          height=10)}) 
      annotation (Placement(transformation(extent={{40,-60},{80,-20}}, rotation=0)));
    Modelica.Blocks.Logical.Hysteresis hysteresis(
      uLow=1.1e5,
      uHigh=2.5e5,
      pre_y_start=true) "mass flow rate signal by pressure control" 
      annotation (Placement(transformation(extent={{-140,-30},{-120,-10}},
            rotation=0)));
    Modelica.Blocks.Logical.Switch switch1 annotation (Placement(transformation(
            extent={{-100,-30},{-80,-10}}, rotation=0)));
    Modelica.Blocks.Sources.Constant m_flow_off(k=0) 
      annotation (Placement(transformation(extent={{-140,10},{-120,30}}, rotation=
             0)));
    Modelica.Blocks.Sources.Constant m_flow_on(k=2) 
      annotation (Placement(transformation(extent={{-140,-60},{-120,-40}},
            rotation=0)));
    Modelica_Fluid.Pipes.StaticPipe overflow(
      redeclare package Medium = Modelica.Media.Water.StandardWater,
      diameter=0.02,
      length=200,
      height_ab=-20) 
                  annotation (Placement(transformation(
          origin={20,-10},
          extent={{10,-10},{-10,10}},
          rotation=90)));
  equation
    connect(massFlowRate.ports[1], upperTank.ports[1]) 
                                                   annotation (Line(
        points={{-40,-30},{-25.3333,-30},{-25.3333,20}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pressure.p, hysteresis.u) annotation (Line(
        points={{61,26},{70,26},{70,80},{-150,80},{-150,-20},{-142,-20}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(hysteresis.y, switch1.u2) annotation (Line(
        points={{-119,-20},{-102,-20}},
        color={255,0,255},
        smooth=Smooth.None));
    connect(m_flow_off.y, switch1.u1) annotation (Line(
        points={{-119,20},{-119,5},{-102,5},{-102,-12}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(m_flow_on.y, switch1.u3) annotation (Line(
        points={{-119,-50},{-110,-50},{-110,-28},{-102,-28}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(switch1.y, massFlowRate.m_flow_in) annotation (Line(
        points={{-79,-20},{-70,-20},{-70,-22},{-60,-22}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(upperTank.ports[2], pipe.port_a) annotation (Line(
        points={{-20,20},{-20,10},{0,10},{0,-20},{6.12323e-016,-20}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_a, pressure.port) annotation (Line(
        points={{6.12323e-016,-20},{0,-20},{0,10},{50,10},{50,16}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_b, lowerTank.ports[1]) annotation (Line(
        points={{-6.12323e-016,-40},{0,-40},{0,-70},{56,-70},{56,-60}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(upperTank.ports[3], overflow.port_a) annotation (Line(
        points={{-14.6667,20},{0,20},{0,40},{20,40},{20,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(overflow.port_b, lowerTank.ports[2]) annotation (Line(
        points={{20,-20},{20,-40},{40,-40},{40,-60},{64,-60}},
        color={0,127,255},
        smooth=Smooth.None));

    annotation (Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-160,-120},{100,120}},
          initialScale=0.1), graphics),
      experiment(StopTime=25000, NumberOfIntervals=5000),
      experimentSetupOutput(equdistant=false),
      Commands(file=
            "../Scripts/Examples/TanksWithOverflow/plot level and port.m_flow.mos"
          "plot level and port.m_flow"),
      Documentation(info="<html>
<p align=justify>The mass flow rate to the upper tank is controlled by the static pressure at its bottom. 
The fluid flows through a pipe and forced by different heights from the upper tank to the lower tank.
</p>
<p>
Additional fluid flows through an overflow pipe if the level of the upper tank exceeds 10m. 
Initially the overflow enters the lower tank above its fluid level; later on the fluid level exceeds the overflow port.
</p>
<p>
Note that the number of solver intervals has been increased, accounting for the long simulation time horizon.
Otherwise the simulation may fail due to too large steps subject to events. Alternatively the 
simulation accuracy could be increased in order to avoid errors. 
</p>
</html>"));
  end TanksWithOverflow;

  model EmptyTanks "Show the treatment of empty tanks"
    extends Modelica.Icons.Example;
    Modelica_Fluid.Vessels.SimpleTank tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      nPorts=1,
      crossArea=1,
      level_start=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=
          0.1)},
      height=1.1)             annotation (Placement(transformation(extent={{-40,20},
              {0,60}},      rotation=0)));

    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=-1) 
                   annotation (Placement(transformation(
          origin={-20,-20},
          extent={{-10,-10},{10,10}},
          rotation=270)));

    Vessels.SimpleTank tank2(
      crossArea=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      nPorts=1,
      height=1.1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.VesselPortsData(diameter=
          0.1, height=0.5)},
      level_start=0) 
      annotation (Placement(transformation(extent={{0,-80},{40,-40}},
            rotation=0)));
    inner Modelica_Fluid.System system 
                                     annotation (Placement(transformation(
            extent={{60,60},{80,80}}, rotation=0)));
  equation
    connect(tank1.ports[1], pipe.port_a) annotation (Line(
        points={{-20,20},{-20,5},{-20,-10},{-20,-10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(pipe.port_b, tank2.ports[1]) annotation (Line(
        points={{-20,-30},{-20,-60},{0,-60},{0,-80},{20,-80}},
        color={0,127,255},
        smooth=Smooth.None));

    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=50),
      experimentSetupOutput,
      Commands(file="../Scripts/Examples/EmptyTanks/plot level and port.p.mos"
          "plot level and port.p"));
  end EmptyTanks;

end Tanks;
