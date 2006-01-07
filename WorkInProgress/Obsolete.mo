package Obsolete "Components that are obsolete and will be removed" 
  extends Modelica.Icons.Library;
  package Types 
    package InitTypes 
      "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
      extends Modelica.Icons.Enumeration;
      constant Integer NoInit = 0 "No explicit initial conditions";
      constant Integer InitialValues = 1 
        "Initial conditions specified by initial values";
      constant Integer SteadyState = 2 "Full steady-state";
      constant Integer SteadyStateHydraulic = 3 
        "Hydraulic steady state, initial values of other state variables given";
      type Temp 
        "Temporary type with choices for menus (until enumerations are available)" 
        extends Integer(min=0,max=3);
        annotation (Evaluate=true, choices(
            choice=Modelica_Fluid.Types.Init.NoInit 
              "NoInit (No explicit initial conditions)",
            choice=Modelica_Fluid.Types.Init.InitialValues 
              "InitialValues (Initial conditions specified by initial values)",
            choice=Modelica_Fluid.Types.Init.SteadyState 
              "SteadyState (Full steady state)",
            choice=Modelica_Fluid.Types.Init.SteadyStateHydraulic 
              "SteadyStateHydraulic (Hydraulic steady state, initial values of other states given)"));
      end Temp;
    end InitTypes;
  end Types;
  
  package Examples 
    model TestPressureLoss 
      extends Modelica.Icons.Example;
      replaceable package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater 
        extends Modelica.Media.Interfaces.PartialMedium 
        "Medium in all components"                        annotation (
        choicesAllMatching =                                                                            true);
      
      parameter WorkInProgress.Utilities.PressureLossFactors lossFactors1(
          zeta1=0.5,
          D_a=0.1,
        zeta2=1) "Loss factor without data for laminar region"   annotation (extent=[-100,-40; -80,-20]);
      parameter WorkInProgress.Utilities.PressureLossFactors lossFactors2(
          zeta1=0.5,
          zetaLaminarKnown=true,
          c0=200,
          D_a=0.1,
        zeta2=1) "Same as lossFactors1 but with data for laminar region" annotation (extent=[-100,-82; -80,-62]);
      
      Sources.FixedAmbient_pTX ambient_p1(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,60; 20,80]);
      Components.PressureLoss orifice_p1(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1) 
        annotation (extent=[-20,60; 0,80]);
      
      annotation (Diagram,
        experiment(StopTime=10, NumberOfIntervals=50000),
        experimentSetupOutput,
        Coordsys(extent=[-100,-200; 200,100]),
        Documentation(info="<html>
<p>
Test whether the different settings of \"from_dp\" and \"use_Re\"
gives the same results for an orifice where the laminar
region is defined and where it is not defined.
</p>
</html>"));
      Sources.PrescribedAmbient_pTX ambient(redeclare package Medium = Medium) 
        annotation (extent=[-60,60; -40,80]);
      Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.999e5; 10,1.001e5]) 
        annotation (extent=[-100,60; -80,80]);
      Sources.FixedAmbient_pTX ambient_p2(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,30; 20,50]);
      Sources.PrescribedMassFlowRate_TX pump_m1(redeclare package Medium = Medium) 
        annotation (extent=[100,60; 120,80]);
      Modelica.Blocks.Sources.TimeTable m_flow_table(table=[0,-10; 10,10]) 
        annotation (extent=[60,60; 80,80]);
      Components.PressureLoss orifice_m1(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1) 
        annotation (extent=[140,60; 160,80]);
      Sources.FixedAmbient_pTX ambient_m1(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,60; 180,80]);
      Components.PressureLoss orifice_p2(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1,
        from_dp=false) 
        annotation (extent=[-20,30; 0,50]);
      
      Sources.FixedAmbient_pTX ambient_p3(
                                         redeclare package Medium = Medium, p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,0; 20,20]);
      Components.PressureLoss orifice_p3(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1,
        use_Re=false) 
        annotation (extent=[-20,0; 0,20]);
      Sources.FixedAmbient_pTX ambient_p4(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-30; 20,-10]);
      Components.PressureLoss orifice_p4(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1,
        from_dp=false,
        use_Re=false) 
        annotation (extent=[-20,-30; 0,-10]);
      Sources.FixedAmbient_pTX ambient_p5(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-60; 20,-40]);
      Components.PressureLoss orifice_p5(
        redeclare package Medium = Medium,
        lossFactors=lossFactors2) 
        annotation (extent=[-20,-60; 0,-40]);
      Sources.FixedAmbient_pTX ambient_p6(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-90; 20,-70]);
      Components.PressureLoss orifice_p6(
        redeclare package Medium = Medium,
        from_dp=false,
        lossFactors=lossFactors2) 
        annotation (extent=[-20,-90; 0,-70]);
      Sources.FixedAmbient_pTX ambient_p7(
                                         redeclare package Medium = Medium, p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-120; 20,-100]);
      Components.PressureLoss orifice_p7(
        redeclare package Medium = Medium,
        use_Re=false,
        lossFactors=lossFactors2) 
        annotation (extent=[-20,-120; 0,-100]);
      Sources.FixedAmbient_pTX ambient_p8(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-150; 20,-130]);
      Components.PressureLoss orifice_p8(
        redeclare package Medium = Medium,
        from_dp=false,
        use_Re=false,
        lossFactors=lossFactors2) 
        annotation (extent=[-20,-150; 0,-130]);
      Sources.PrescribedMassFlowRate_TX pump_m2(redeclare package Medium = Medium) 
        annotation (extent=[100,30; 120,50]);
      Components.PressureLoss orifice_m2(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1,
        from_dp=false) 
        annotation (extent=[140,30; 160,50]);
      Sources.FixedAmbient_pTX ambient_m2(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,30; 180,50]);
      Sources.PrescribedMassFlowRate_TX pump_m3(redeclare package Medium = Medium) 
        annotation (extent=[100,0; 120,20]);
      Components.PressureLoss orifice_m3(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1,
        use_Re=false) 
        annotation (extent=[140,0; 160,20]);
      Sources.FixedAmbient_pTX ambient_m3(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,0; 180,20]);
      Sources.PrescribedMassFlowRate_TX pump_m4(redeclare package Medium = Medium) 
        annotation (extent=[100,-30; 120,-10]);
      Components.PressureLoss orifice_m4(
        redeclare package Medium = Medium,
        lossFactors=lossFactors1,
        from_dp=false,
        use_Re=false) 
        annotation (extent=[140,-30; 160,-10]);
      Sources.FixedAmbient_pTX ambient_m4(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-30; 180,-10]);
      Sources.PrescribedMassFlowRate_TX pump_m5(redeclare package Medium = Medium) 
        annotation (extent=[100,-60; 120,-40]);
      Components.PressureLoss orifice_m5(
        redeclare package Medium = Medium,
        lossFactors=lossFactors2) 
        annotation (extent=[140,-60; 160,-40]);
      Sources.FixedAmbient_pTX ambient_m5(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-60; 180,-40]);
      Sources.PrescribedMassFlowRate_TX pump_m6(redeclare package Medium = Medium) 
        annotation (extent=[100,-90; 120,-70]);
      Components.PressureLoss orifice_m6(
        redeclare package Medium = Medium,
        from_dp=false,
        lossFactors=lossFactors2) 
        annotation (extent=[140,-90; 160,-70]);
      Sources.FixedAmbient_pTX ambient_m6(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-90; 180,-70]);
      Sources.PrescribedMassFlowRate_TX pump_m7(redeclare package Medium = Medium) 
        annotation (extent=[100,-120; 120,-100]);
      Components.PressureLoss orifice_m7(
        redeclare package Medium = Medium,
        use_Re=false,
        lossFactors=lossFactors2) 
        annotation (extent=[140,-120; 160,-100]);
      Sources.FixedAmbient_pTX ambient_m7(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-120; 180,-100]);
      Sources.PrescribedMassFlowRate_TX pump_m8(redeclare package Medium = Medium) 
        annotation (extent=[100,-150; 120,-130]);
      Components.PressureLoss orifice_m8(
        redeclare package Medium = Medium,
        from_dp=false,
        use_Re=false,
        lossFactors=lossFactors2) 
        annotation (extent=[140,-150; 160,-130]);
      Sources.FixedAmbient_pTX ambient_m8(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-150; 180,-130]);
    equation 
      connect(orifice_p1.port_b, ambient_p1.port) 
        annotation (points=[1,70; 20,70], style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, orifice_p1.port_a) 
                                             annotation (points=[-40,70; -21,
            70],
          style(color=69, rgbcolor={0,127,255}));
      connect(p_table.y, ambient.p_in)  annotation (points=[-79,70; -72,70; -72,76;
            -62,76], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m1.m_flow_in) annotation (points=[81,70; 88,70;
            88,76; 100.7,76], style(color=74, rgbcolor={0,0,127}));
      connect(pump_m1.port, orifice_m1.port_a) annotation (points=[120,70;
            139,70],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m1.port_b, ambient_m1.port) annotation (points=[161,70;
            180,70],
                 style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, orifice_p2.port_a) annotation (points=[-40,70;
            -30,70; -30,40; -21,40],
                             style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p2.port_b, ambient_p2.port) 
        annotation (points=[1,40; 20,40], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p3.port_b, ambient_p3.port) 
        annotation (points=[1,10; 20,10], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p4.port_b, ambient_p4.port) 
        annotation (points=[1,-20; 20,-20], style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, orifice_p3.port_a) annotation (points=[-40,70;
            -30,70; -30,10; -21,10],
                             style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, orifice_p4.port_a) annotation (points=[-40,70;
            -30,70; -30,-20; -21,-20],
                               style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p5.port_b, ambient_p5.port) 
        annotation (points=[1,-50; 20,-50],
                                          style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port,orifice_p5. port_a) 
                                             annotation (points=[-40,70; -30,
            70; -30,-50; -21,-50],
          style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, orifice_p6.port_a) annotation (points=[-40,70;
            -30,70; -30,-80; -21,-80],
                               style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p6.port_b, ambient_p6.port) 
        annotation (points=[1,-80; 20,-80], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p7.port_b, ambient_p7.port) 
        annotation (points=[1,-110; 20,-110],
                                          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p8.port_b, ambient_p8.port) annotation (points=[1,-140;
            20,-140],
                   style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, orifice_p7.port_a) annotation (points=[-40,70;
            -30,70; -30,-110; -21,-110],
                                 style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, orifice_p8.port_a) annotation (points=[-40,70;
            -30,70; -30,-140; -21,-140],
                                 style(color=69, rgbcolor={0,127,255}));
      connect(pump_m2.port, orifice_m2.port_a) annotation (points=[120,40;
            139,40],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m2.port_b, ambient_m2.port) annotation (points=[161,40;
            180,40],
                 style(color=69, rgbcolor={0,127,255}));
      connect(pump_m3.port, orifice_m3.port_a) annotation (points=[120,10;
            139,10],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m3.port_b, ambient_m3.port) annotation (points=[161,10;
            180,10],
                 style(color=69, rgbcolor={0,127,255}));
      connect(pump_m4.port, orifice_m4.port_a) annotation (points=[120,-20;
            139,-20],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m4.port_b, ambient_m4.port) annotation (points=[161,-20;
            180,-20],
                  style(color=69, rgbcolor={0,127,255}));
      connect(m_flow_table.y, pump_m2.m_flow_in) annotation (points=[81,70; 88,70;
            88,46; 100.7,46], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m3.m_flow_in) annotation (points=[81,70; 88,70;
            88,16; 100.7,16], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m4.m_flow_in) annotation (points=[81,70; 88,70;
            88,-14; 100.7,-14], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m5.m_flow_in) annotation (points=[81,70; 88,70;
            88,-44; 100.7,-44], style(color=74, rgbcolor={0,0,127}));
      connect(pump_m5.port, orifice_m5.port_a) annotation (points=[120,-50;
            139,-50],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m5.port_b, ambient_m5.port) annotation (points=[161,-50;
            180,-50],
                  style(color=69, rgbcolor={0,127,255}));
      connect(pump_m6.port, orifice_m6.port_a) annotation (points=[120,-80;
            139,-80],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m6.port_b, ambient_m6.port) annotation (points=[161,-80;
            180,-80],
                  style(color=69, rgbcolor={0,127,255}));
      connect(pump_m7.port, orifice_m7.port_a) annotation (points=[120,-110;
            139,-110],
                   style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m7.port_b, ambient_m7.port) annotation (points=[161,-110;
            180,-110],
                   style(color=69, rgbcolor={0,127,255}));
      connect(pump_m8.port, orifice_m8.port_a) annotation (points=[120,-140;
            139,-140],
                   style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m8.port_b, ambient_m8.port) annotation (points=[161,-140;
            180,-140],
                   style(color=69, rgbcolor={0,127,255}));
      connect(m_flow_table.y, pump_m6.m_flow_in) annotation (points=[81,70; 88,70;
            88,-74; 100.7,-74], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m7.m_flow_in) annotation (points=[81,70; 88,70;
            88,-104; 100.7,-104], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m8.m_flow_in) annotation (points=[81,70; 88,70;
            88,-134; 100.7,-134], style(color=74, rgbcolor={0,0,127}));
    end TestPressureLoss;
    
    model TestSharpEdgedOrifice 
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      // Modelica.Media.Water.ConstantPropertyLiquidWater 
      // Modelica.Media.Water.WaterIF97_ph 
      replaceable package Medium = Modelica.Media.IdealGases.SingleGases.O2 
        extends Modelica.Media.Interfaces.PartialMedium 
        "Medium in all components"                        annotation (
        choicesAllMatching =                                                                            true);
      
      annotation (Diagram,
        experiment(StopTime=10, NumberOfIntervals=10000),
        experimentSetupOutput,
        Coordsys(extent=[-100,-100; 100,100]),
        Documentation(info="<html>
</html>"));
      Sources.PrescribedAmbient_pTX ambient(redeclare package Medium = Medium) 
        annotation (extent=[-40,40; -20,60]);
      Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.1e5; 10,2e5]) 
        annotation (extent=[-80,40; -60,60]);
      
      Sources.FixedAmbient_pTX ambient_p1(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[60,40; 40,60]);
      Sources.FixedAmbient_pTX ambient_p2(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[60,10; 40,30]);
      Modelica_Fluid.WorkInProgress.Components.SharpEdgedOrifice orifice1(
        redeclare package Medium = Medium,
        D_pipe=0.1,
        alpha=50,
        D_min=0.02,
        L=0.005) annotation (extent=[0,40; 20,60]);
      Modelica_Fluid.WorkInProgress.Components.SharpEdgedOrifice orifice2(
        redeclare package Medium = Medium,
        D_pipe=0.1,
        alpha=50,
        from_dp=false,
        D_min=0.02,
        L=0.005) annotation (extent=[0,10; 20,30]);
    equation 
      connect(p_table.y, ambient.p_in)  annotation (points=[-59,50; -52,50; -52,56;
            -42,56], style(color=74, rgbcolor={0,0,127}));
      connect(ambient.port, orifice1.port_a) annotation (points=[-20,50; -1,
            50],
          style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8));
      connect(orifice1.port_b, ambient_p1.port) annotation (points=[21,50; 40,
            50],
          style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8));
      connect(ambient.port, orifice2.port_a) annotation (points=[-20,50; -12,
            50; -12,20; -1,20],
                            style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8));
      connect(orifice2.port_b, ambient_p2.port) annotation (points=[21,20; 40,
            20],
          style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=3,
          rgbfillColor={0,0,255},
          fillPattern=8));
    end TestSharpEdgedOrifice;
    
    model TestSuddenExpansion 
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      replaceable package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater 
        extends Modelica.Media.Interfaces.PartialMedium 
        "Medium in all components"                        annotation (
        choicesAllMatching =                                                                            true);
      
      annotation (Diagram,
        experiment(StopTime=10, NumberOfIntervals=10000),
        experimentSetupOutput,
        Coordsys(extent=[-100,-100; 100,100]),
        Documentation(info="<html>
</html>"));
      Sources.PrescribedAmbient_pTX ambient(redeclare package Medium = Medium) 
        annotation (extent=[-40,40; -20,60]);
      Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.9999e5; 10,1.0001e5]) 
        annotation (extent=[-80,40; -60,60]);
      
      Sources.FixedAmbient_pTX ambient_p1(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[62,40; 42,60]);
      Sources.FixedAmbient_pTX ambient_p2(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[60,10; 40,30]);
      Modelica_Fluid.WorkInProgress.Components.SuddenExpansion expansion1(
        redeclare package Medium = Medium,
        D_a=0.1,
        D_b=0.2) annotation (extent=[0,40; 20,60]);
      Modelica_Fluid.WorkInProgress.Components.SuddenExpansion expansion2(
        redeclare package Medium = Medium,
        D_a=0.1,
        D_b=0.2,
        from_dp=false) annotation (extent=[0,10; 20,30]);
    equation 
      connect(p_table.y, ambient.p_in)  annotation (points=[-59,50; -52,50; -52,56;
            -42,56], style(color=74, rgbcolor={0,0,127}));
      connect(ambient.port, expansion1.port_a) annotation (points=[-20,50; -1,
            50],
          style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(expansion1.port_b, ambient_p1.port) annotation (points=[21,50;
            42,50],
          style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(expansion2.port_b, ambient_p2.port) annotation (points=[21,20;
            40,20],
          style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
      connect(expansion2.port_a, ambient.port) annotation (points=[-1,20; -10,
            20; -10,50; -20,50],
                             style(
          color=69,
          rgbcolor={0,127,255},
          fillColor=7,
          rgbfillColor={255,255,255},
          fillPattern=1));
    end TestSuddenExpansion;
    
    model TestWallFriction 
      import Modelica_Fluid;
      extends Modelica.Icons.Example;
      replaceable package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater 
        extends Modelica.Media.Interfaces.PartialMedium 
        "Medium in all components"                        annotation (
        choicesAllMatching =                                                                            true);
      parameter Modelica.SIunits.Length roughness = 0.025e-3;
      
      annotation (Diagram,
        experiment(StopTime=10, NumberOfIntervals=10000),
        experimentSetupOutput,
        Coordsys(extent=[-100,-100; 100,100]),
        Documentation(info="<html>
</html>"));
      Sources.PrescribedAmbient_pTX ambient(redeclare package Medium = Medium) 
        annotation (extent=[-40,40; -20,60]);
      Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.9999e5; 10,1.0001e5]) 
        annotation (extent=[-80,40; -60,60]);
      
      Sources.FixedAmbient_pTX ambient_p1(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[62,40; 42,60]);
      Sources.FixedAmbient_pTX ambient_p2(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[60,10; 40,30]);
      Components.WallFriction pipe1(
        length=1,
        diameter=0.1,
        redeclare package Medium = Medium,
        roughness=roughness) annotation (extent=[0,40; 20,60]);
      Components.WallFriction pipe2(
        length=1,
        diameter=0.1,
        from_dp=false,
        redeclare package Medium = Medium,
        roughness=roughness) annotation (extent=[0,10; 20,30]);
      Modelica_Fluid.Components.PressureDropPipe pressureDropPipe(
        redeclare package Medium = Medium,
        frictionType=Modelica_Fluid.Types.FrictionTypes.DetailedFriction,
        length=1,
        diameter=0.1,
        roughness=roughness) annotation (extent=[0,-20; 20,0]);
      Sources.FixedAmbient_pTX ambient_p3(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[60,-20; 40,0]);
    equation 
      connect(p_table.y, ambient.p_in)  annotation (points=[-59,50; -52,50; -52,56;
            -42,56], style(color=74, rgbcolor={0,0,127}));
      connect(ambient.port, pipe1.port_a) annotation (points=[-20,50; -1,50],
          style(color=69, rgbcolor={0,127,255}));
      connect(pipe1.port_b, ambient_p1.port) 
        annotation (points=[21,50; 42,50],style(color=69, rgbcolor={0,127,255}));
      connect(pipe2.port_b, ambient_p2.port) 
        annotation (points=[21,20; 40,20],style(color=69, rgbcolor={0,127,255}));
      connect(ambient.port, pipe2.port_a) annotation (points=[-20,50; -12,50; -12,
            20; -1,20],  style(color=69, rgbcolor={0,127,255}));
      connect(pressureDropPipe.port_b, ambient_p3.port) 
        annotation (points=[20,-10; 40,-10],
                                          style(color=69, rgbcolor={0,127,255}));
      connect(pressureDropPipe.port_a, ambient.port)  annotation (points=[0,-10;
            -12,-10; -12,50; -20,50],style(color=69, rgbcolor={0,127,255}));
    end TestWallFriction;
  end Examples;
end Obsolete;
