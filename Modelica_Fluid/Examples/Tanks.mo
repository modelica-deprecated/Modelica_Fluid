within Modelica_Fluid.Examples;
package Tanks "Library demonstrating the usage of the tank model"
  extends Modelica.Icons.Library;
  model OneTank
    "Demonstrates a tank with one constant top inlet mass flow rate and a bottom outlet into the ambient"
    import Modelica.SIunits.Conversions.from_bar;
    extends Modelica.Icons.Example;

    Modelica_Fluid.Vessels.Tank tank(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      height=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)},
      level_start=0.9,
      V0=0.1,
      nTopPorts=1,
      nPorts=1,
      stiffCharacteristicForEmptyPort=true) 
      annotation (Placement(transformation(extent={{-40,28},{0,68}}, rotation=0)));

    Sources.MassFlowSource_T flowSource(nPorts=1,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      m_flow=20,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-52,70},{-32,90}}, rotation=
             0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=100),
      experimentSetupOutput,
      Commands(file="../Scripts/Examples/OneTank/plot level and port.p.mos"
          "plot level and port.p", file=
            "../Scripts/Examples/OneTank/plot level, port.p and port.m_flow.mos"
          "plot level, port.p and port.m_flow"));
    inner Modelica_Fluid.System system 
                          annotation (Placement(transformation(extent={{-10,72},
              {10,92}}, rotation=0)));
    Modelica_Fluid.Sources.Boundary_pT ambient_fixed(nPorts=1,
                                           redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      p=system.p_ambient,
      T=system.T_ambient) 
      annotation (Placement(transformation(extent={{-54,-20},{-34,0}}, rotation=
             0)));
    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={-20,10},
          extent={{-10,-10},{10,10}},
          rotation=90)));
  equation
    connect(flowSource.ports[1], tank.topPorts[1])  annotation (Line(points={{-32,80},
            {-19,80},{-19,68}},     color={0,127,255}));
    connect(ambient_fixed.ports[1], pipe.port_a) annotation (Line(points={{-34,-10},
            {-20,-10},{-20,0}}, color={0,127,255}));
    connect(pipe.port_b, tank.ports[1]) annotation (Line(points={{-20,20},{-20,
            24},{-20,28},{-21,28}},
                  color={0,127,255}));
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
    Modelica_Fluid.Vessels.Tank tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort,
      crossArea=1,
      height=4,
      level_start=3,
      T_start=Modelica.SIunits.Conversions.from_degC(50),
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)}) 
      annotation (Placement(transformation(extent={{-80,0},{-40,40}}, rotation=
              0)));
    Modelica_Fluid.Vessels.Tank tank2(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort,
      crossArea=1,
      height=4,
      level_start=1,
      T_start=Modelica.SIunits.Conversions.from_degC(100),
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)}) 
      annotation (Placement(transformation(extent={{0,0},{40,40}}, rotation=0)));
    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1)  annotation (Placement(transformation(extent={{-30,-30},{
              -10,-10}}, rotation=0)));
  equation
    connect(tank1.ports[1], pipe.port_a) annotation (Line(points={{-61,0},{-61,
            -20},{-30,-20}}, color={0,127,255}));
    connect(pipe.port_b, tank2.ports[1]) annotation (Line(points={{-10,-20},{19,
            -20},{19,0}},  color={0,127,255}));
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
    Modelica_Fluid.Vessels.Tank tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      level_start=0.1,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.1,
          portLevel=1)},
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
            {-20,-30},{-20,-22},{-21,-22}},
                             color={0,127,255}));
    connect(pipe.port_a, tank1.ports[2]) annotation (Line(points={{40,0},{40,
            -28},{-18,-28},{-18,-20},{-21,-20},{-21,-18}}, color={0,127,255}));
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
    Modelica_Fluid.Vessels.Tank tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.1,
          portLevel=1)},
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
    connect(tank1.ports[1], pipe1.port_b) annotation (Line(points={{-21,-22},{
            -21,-35},{-20,-35},{-20,-50}},
                       color={0,127,255}));
    connect(ambient_fixed.ports[1], pipe1.port_a) annotation (Line(points={{-40,-90},
            {-20,-90},{-20,-70}}, color={0,127,255}));
    connect(tank1.ports[2], pipe2.port_b) annotation (Line(points={{-21,-18},{
            -18,-18},{-18,-40},{30,-40},{30,-50}}, color={0,127,255}));
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
    Modelica_Fluid.Vessels.Tank tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.1,
          portLevel=1)},
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
    Modelica_Fluid.Vessels.Tank tank2(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      crossArea=1,
      V0=0.1,
      height=2,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.1,
          portLevel=0.5)},
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
    connect(tank1.ports[1], pipe1.port_b) annotation (Line(points={{-61,-2},{
            -61,-15},{-60,-15},{-60,-30}},
                       color={0,127,255}));
    connect(ambient_fixed1.ports[1], pipe1.port_a) 
                                              annotation (Line(points={{-80,-70},
            {-60,-70},{-60,-50}}, color={0,127,255}));
    connect(ambient_fixed2.ports[1], pipe2.port_a) annotation (Line(points={{20,-70},
            {40,-70},{40,-50}}, color={0,127,255}));
    connect(tank2.ports[1], pipe2.port_b) 
      annotation (Line(points={{39,-2},{39,-15},{40,-15},{40,-30}},
                                                  color={0,127,255}));
    connect(pipe3.port_a, tank1.ports[2]) annotation (Line(points={{-20,20},{
            -30,20},{-30,-10},{-58,-10},{-58,0},{-61,0},{-61,2}},  color={0,127,
            255}));
    connect(pipe3.port_b, tank2.ports[2]) annotation (Line(points={{0,20},{10,
            20},{10,-8},{38,-8},{38,0},{39,0},{39,2}},  color={0,127,255}));
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
    Modelica_Fluid.Vessels.Tank tank3(
      redeclare package Medium = Medium,
      crossArea=1,
      V0=0.1,
      height=20,
      nPorts=2,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.05,
          portLevel=6.5)},
      level_start=6,
      nTopPorts=1,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{-80,-50},{-40,-10}},
            rotation=0)));
    Modelica_Fluid.Vessels.Tank tank1(
      redeclare package Medium = Medium,
      crossArea=1,
      V0=0.1,
      height=10,
      nPorts=1,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.1, portLevel=0)},
      level_start=9,
      stiffCharacteristicForEmptyPort = stiffCharacteristicForEmptyPort) 
      annotation (Placement(transformation(extent={{50,50},{90,90}}, rotation=0)));
    Modelica_Fluid.Vessels.Tank tank2(
      redeclare package Medium = Medium,
      crossArea=1,
      V0=0.1,
      height=10,
      nPorts=3,
      portsData={Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.05, portLevel=0),
          Modelica_Fluid.Vessels.BaseClasses.TankPortData(diameter=0.05,
          portLevel=2),Modelica_Fluid.Vessels.BaseClasses.TankPortData(
          diameter=0.1, portLevel=3)},
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
            {-60,-59},{-60,-52},{-61,-52}},
                             color={0,127,255}));
    connect(pipe1.port_b, tank1.ports[1]) annotation (Line(points={{70,40},{70,
            45},{70,50},{69,50}},
                  color={0,127,255}));
    connect(pipe2.port_a, tank3.ports[2]) annotation (Line(points={{
            -6.12323e-016,-32},{-6.12323e-016,-48},{0,-60},{-58,-60},{-58,-48},
            {-61,-48}}, color={0,127,255}));
    connect(pipe3.port_a, tank3.topPorts[1]) 
                                            annotation (Line(points={{-60,0},{
            -60,-5},{-60,-10},{-59,-10}},
                      color={0,127,255}));
    connect(pipe3.port_b, tank2.ports[1]) annotation (Line(points={{-60,20},{
            -60,26},{-30,26},{-30,0},{-2,0},{-2,7.33333},{-1,7.33333}},
                                                           color={0,127,255}));
    connect(pipe1.port_a, tank2.ports[3]) annotation (Line(points={{70,20},{70,
            0},{2,0},{2,12.6667},{-1,12.6667}},
                                   color={0,127,255}));
    connect(pipe2.port_b, tank2.ports[2]) annotation (Line(
        points={{6.12323e-016,-12},{0,-12},{0,10},{-1,10}},
        color={0,127,255},
        smooth=Smooth.None));
  end TanksWithEmptyingPipe2;

  model ThreeOpenTanks "Demonstrating the usage of OpenTank"
    import Modelica_Fluid;
    extends Modelica.Icons.Example;
     // replaceable package Medium = Modelica_Fluid.Media.Water.ConstantPropertyLiquidWater extends
    // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends
    // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
     replaceable package Medium = 
        Modelica.Media.Water.ConstantPropertyLiquidWater                           constrainedby
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    Modelica_Fluid.Vessels.OpenTank tank1(
      crossArea=1,
      redeclare package Medium = Medium,
      use_portDiameters=true,
      height=12,
      level_start=8,
      nPorts=1,
      portDiameters={0.1}) 
                     annotation (Placement(transformation(extent={{-80,20},{-40,
              60}}, rotation=0)));
    Modelica_Fluid.Vessels.OpenTank tank2(
      crossArea=1,
      redeclare package Medium = Medium,
      use_portDiameters=true,
      height=12,
      level_start=3,
      nPorts=1,
      portDiameters={0.1}) 
                     annotation (Placement(transformation(extent={{-20,20},{20,
              60}}, rotation=0)));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
              -100},{100,100}}),
                        graphics),
      experiment(StopTime=200),
      experimentSetupOutput,
      Documentation(info="<html>
  
</html>"));

    inner Modelica_Fluid.System system(energyDynamics=Modelica_Fluid.Types.Dynamics.FixedInitial) 
                                     annotation (Placement(transformation(
            extent={{70,-90},{90,-70}}, rotation=0)));
    Modelica_Fluid.Vessels.OpenTank tank3(
      crossArea=1,
      redeclare package Medium = Medium,
      use_portDiameters=true,
      height=12,
      level_start=3,
      nPorts=1,
      portDiameters={0.1}) 
                     annotation (Placement(transformation(extent={{40,10},{80,50}},
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
  end ThreeOpenTanks;

  model TestEmptyOpenTank "Test whether an empty tank is properly handeled"
    extends Modelica.Icons.Example;
    Vessels.OpenTank tank1(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      nPorts=1,
      crossArea=1,
      level_start=1,
      height=1.1,
      portDiameters={0.1},
      V0=1e-3)                annotation (Placement(transformation(extent={{-20,
              20},{20,60}}, rotation=0)));

    Modelica_Fluid.Pipes.StaticPipe pipe(
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      length=1,
      diameter=0.1,
      height_ab=1) annotation (Placement(transformation(
          origin={0,-10},
          extent={{10,-10},{-10,10}},
          rotation=270)));

    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      experiment(StopTime=50),
      experimentSetupOutput);
    Modelica_Fluid.Vessels.Tank tank2(
      nTopPorts=1,
      nPorts=0,
      crossArea=1,
      level_start=0,
      redeclare package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater,
      height=1.1,
      V0=1e-3) 
      annotation (Placement(transformation(extent={{-20,-80},{20,-40}},
            rotation=0)));
    inner Modelica_Fluid.System system 
                                     annotation (Placement(transformation(
            extent={{56,58},{76,78}}, rotation=0)));
  equation
    connect(pipe.port_b, tank1.ports[1]) annotation (Line(points={{1.83697e-015,
            0},{0,0},{0,20}}, color={0,127,255}));
    connect(pipe.port_a, tank2.topPorts[1]) annotation (Line(points={{
            -1.83697e-015,-20},{0,-20},{0,-39}}, color={0,127,255}));
  end TestEmptyOpenTank;

end Tanks;
