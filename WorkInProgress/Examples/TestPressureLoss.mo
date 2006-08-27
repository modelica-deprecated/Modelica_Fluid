model TestPressureLoss 
  extends Modelica.Icons.Example;
  replaceable package Medium = 
      Modelica.Media.Water.ConstantPropertyLiquidWater 
    extends Modelica.Media.Interfaces.PartialMedium "Medium in all components" 
                                                      annotation (
    choicesAllMatching =                                                                            true);
  
  parameter 
    Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.LossFactorData
    lossFactors1(
      zeta1=0.5,
      D_a=0.1,
      D_b=0.1,
      Re_turbulent=100,
      D_Re=0.1,
    zeta2=1) "Loss factor without data for laminar region"   annotation (extent=[-100,-40; -80,-20]);
  parameter 
    Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.LossFactorData
    lossFactors2(
      zeta1=0.5,
      zetaLaminarKnown=true,
      c0=200,
      D_a=0.1,
      D_b=0.1,
      Re_turbulent=100,
      D_Re=0.1,
    zeta2=1) "Same as lossFactors1 but with data for laminar region" annotation (extent=[-100,-82; -80,-62]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,60; 20,80]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p1(
    redeclare package Medium = Medium,
    data=lossFactors1) 
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
  Modelica_Fluid.Sources.PrescribedBoundary_pTX ambientSource(
                                              redeclare package Medium = Medium) 
    annotation (extent=[-60,60; -40,80]);
  Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.999e5; 10,1.001e5]) 
    annotation (extent=[-100,60; -80,80]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,30; 20,50]);
  Sources.PrescribedMassFlowRate_TX pump_m1(redeclare package Medium = Medium) 
    annotation (extent=[100,60; 120,80]);
  Modelica.Blocks.Sources.TimeTable m_flow_table(table=[0,-10; 10,10]) 
    annotation (extent=[60,60; 80,80]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m1(
    redeclare package Medium = Medium,
    data=lossFactors1) 
    annotation (extent=[140,60; 160,80]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m1(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,60; 180,80]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p2(
    redeclare package Medium = Medium,
    data=lossFactors1,
    from_dp=false) 
    annotation (extent=[-20,30; 0,50]);
  
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p3(
                                     redeclare package Medium = Medium, p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,0; 20,20]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p3(
    redeclare package Medium = Medium,
    data=lossFactors1,
    use_Re=false) 
    annotation (extent=[-20,0; 0,20]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p4(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,-30; 20,-10]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p4(
    redeclare package Medium = Medium,
    data=lossFactors1,
    from_dp=false,
    use_Re=false) 
    annotation (extent=[-20,-30; 0,-10]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p5(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,-60; 20,-40]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p5(
    redeclare package Medium = Medium,
    data=lossFactors2) 
    annotation (extent=[-20,-60; 0,-40]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p6(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,-90; 20,-70]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p6(
    redeclare package Medium = Medium,
    from_dp=false,
    data=lossFactors2) 
    annotation (extent=[-20,-90; 0,-70]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p7(
                                     redeclare package Medium = Medium, p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,-120; 20,-100]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p7(
    redeclare package Medium = Medium,
    use_Re=false,
    data=lossFactors2) 
    annotation (extent=[-20,-120; 0,-100]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p8(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[40,-150; 20,-130]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_p8(
    redeclare package Medium = Medium,
    from_dp=false,
    use_Re=false,
    data=lossFactors2) 
    annotation (extent=[-20,-150; 0,-130]);
  Sources.PrescribedMassFlowRate_TX pump_m2(redeclare package Medium = Medium) 
    annotation (extent=[100,30; 120,50]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m2(
    redeclare package Medium = Medium,
    data=lossFactors1,
    from_dp=false) 
    annotation (extent=[140,30; 160,50]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m2(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,30; 180,50]);
  Sources.PrescribedMassFlowRate_TX pump_m3(redeclare package Medium = Medium) 
    annotation (extent=[100,0; 120,20]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m3(
    redeclare package Medium = Medium,
    data=lossFactors1,
    use_Re=false) 
    annotation (extent=[140,0; 160,20]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m3(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,0; 180,20]);
  Sources.PrescribedMassFlowRate_TX pump_m4(redeclare package Medium = Medium) 
    annotation (extent=[100,-30; 120,-10]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m4(
    redeclare package Medium = Medium,
    data=lossFactors1,
    from_dp=false,
    use_Re=false) 
    annotation (extent=[140,-30; 160,-10]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m4(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,-30; 180,-10]);
  Sources.PrescribedMassFlowRate_TX pump_m5(redeclare package Medium = Medium) 
    annotation (extent=[100,-60; 120,-40]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m5(
    redeclare package Medium = Medium,
    data=lossFactors2) 
    annotation (extent=[140,-60; 160,-40]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m5(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,-60; 180,-40]);
  Sources.PrescribedMassFlowRate_TX pump_m6(redeclare package Medium = Medium) 
    annotation (extent=[100,-90; 120,-70]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m6(
    redeclare package Medium = Medium,
    from_dp=false,
    data=lossFactors2) 
    annotation (extent=[140,-90; 160,-70]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m6(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,-90; 180,-70]);
  Sources.PrescribedMassFlowRate_TX pump_m7(redeclare package Medium = Medium) 
    annotation (extent=[100,-120; 120,-100]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m7(
    redeclare package Medium = Medium,
    use_Re=false,
    data=lossFactors2) 
    annotation (extent=[140,-120; 160,-100]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m7(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,-120; 180,-100]);
  Sources.PrescribedMassFlowRate_TX pump_m8(redeclare package Medium = Medium) 
    annotation (extent=[100,-150; 120,-130]);
  Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
    orifice_m8(
    redeclare package Medium = Medium,
    from_dp=false,
    use_Re=false,
    data=lossFactors2) 
    annotation (extent=[140,-150; 160,-130]);
  Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m8(
    redeclare package Medium = Medium,
    p=1.0e5,
    T=Modelica.SIunits.Conversions.from_degC(80)) 
    annotation (extent=[200,-150; 180,-130]);
  inner Ambient ambient annotation (extent=[-100,-140; -80,-120]);
equation 
  connect(orifice_p1.port_b, ambientSource_p1.port) 
    annotation (points=[0,70; 20,70], style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port, orifice_p1.port_a) 
                                         annotation (points=[-40,70; -20,70],
      style(color=69, rgbcolor={0,127,255}));
  connect(p_table.y, ambientSource.p_in)  annotation (points=[-79,70; -72,70; -72,76;
        -62,76], style(color=74, rgbcolor={0,0,127}));
  connect(m_flow_table.y, pump_m1.m_flow_in) annotation (points=[81,70; 88,70;
        88,76; 100.7,76], style(color=74, rgbcolor={0,0,127}));
  connect(pump_m1.port, orifice_m1.port_a) annotation (points=[120,70; 140,70],
      style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m1.port_b, ambientSource_m1.port) annotation (points=[160,70; 180,
        70], style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port, orifice_p2.port_a) annotation (points=[-40,70; -30,70;
        -30,40; -20,40], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_p2.port_b, ambientSource_p2.port) 
    annotation (points=[0,40; 20,40], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_p3.port_b, ambientSource_p3.port) 
    annotation (points=[0,10; 20,10], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_p4.port_b, ambientSource_p4.port) 
    annotation (points=[0,-20; 20,-20], style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port, orifice_p3.port_a) annotation (points=[-40,70; -30,70;
        -30,10; -20,10], style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port, orifice_p4.port_a) annotation (points=[-40,70; -30,70;
        -30,-20; -20,-20], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_p5.port_b, ambientSource_p5.port) 
    annotation (points=[0,-50; 20,-50],
                                      style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port,orifice_p5. port_a) 
                                         annotation (points=[-40,70; -30,70;
        -30,-50; -20,-50],
      style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port, orifice_p6.port_a) annotation (points=[-40,70; -30,70;
        -30,-80; -20,-80], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_p6.port_b, ambientSource_p6.port) 
    annotation (points=[0,-80; 20,-80], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_p7.port_b, ambientSource_p7.port) 
    annotation (points=[0,-110; 20,-110],
                                      style(color=69, rgbcolor={0,127,255}));
  connect(orifice_p8.port_b, ambientSource_p8.port) annotation (points=[0,-140; 20,
        -140], style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port, orifice_p7.port_a) annotation (points=[-40,70; -30,70;
        -30,-110; -20,-110], style(color=69, rgbcolor={0,127,255}));
  connect(ambientSource.port, orifice_p8.port_a) annotation (points=[-40,70; -30,70;
        -30,-140; -20,-140], style(color=69, rgbcolor={0,127,255}));
  connect(pump_m2.port, orifice_m2.port_a) annotation (points=[120,40; 140,40],
      style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m2.port_b, ambientSource_m2.port) annotation (points=[160,40; 180,
        40], style(color=69, rgbcolor={0,127,255}));
  connect(pump_m3.port, orifice_m3.port_a) annotation (points=[120,10; 140,10],
      style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m3.port_b, ambientSource_m3.port) annotation (points=[160,10; 180,
        10], style(color=69, rgbcolor={0,127,255}));
  connect(pump_m4.port, orifice_m4.port_a) annotation (points=[120,-20; 140,-20],
      style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m4.port_b, ambientSource_m4.port) annotation (points=[160,-20; 180,
        -20], style(color=69, rgbcolor={0,127,255}));
  connect(m_flow_table.y, pump_m2.m_flow_in) annotation (points=[81,70; 88,70;
        88,46; 100.7,46], style(color=74, rgbcolor={0,0,127}));
  connect(m_flow_table.y, pump_m3.m_flow_in) annotation (points=[81,70; 88,70;
        88,16; 100.7,16], style(color=74, rgbcolor={0,0,127}));
  connect(m_flow_table.y, pump_m4.m_flow_in) annotation (points=[81,70; 88,70;
        88,-14; 100.7,-14], style(color=74, rgbcolor={0,0,127}));
  connect(m_flow_table.y, pump_m5.m_flow_in) annotation (points=[81,70; 88,70;
        88,-44; 100.7,-44], style(color=74, rgbcolor={0,0,127}));
  connect(pump_m5.port, orifice_m5.port_a) annotation (points=[120,-50; 140,-50],
      style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m5.port_b, ambientSource_m5.port) annotation (points=[160,-50; 180,
        -50], style(color=69, rgbcolor={0,127,255}));
  connect(pump_m6.port, orifice_m6.port_a) annotation (points=[120,-80; 140,-80],
      style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m6.port_b, ambientSource_m6.port) annotation (points=[160,-80; 180,
        -80], style(color=69, rgbcolor={0,127,255}));
  connect(pump_m7.port, orifice_m7.port_a) annotation (points=[120,-110; 140,
        -110], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m7.port_b, ambientSource_m7.port) annotation (points=[160,-110; 180,
        -110], style(color=69, rgbcolor={0,127,255}));
  connect(pump_m8.port, orifice_m8.port_a) annotation (points=[120,-140; 140,
        -140], style(color=69, rgbcolor={0,127,255}));
  connect(orifice_m8.port_b, ambientSource_m8.port) annotation (points=[160,-140; 180,
        -140], style(color=69, rgbcolor={0,127,255}));
  connect(m_flow_table.y, pump_m6.m_flow_in) annotation (points=[81,70; 88,70;
        88,-74; 100.7,-74], style(color=74, rgbcolor={0,0,127}));
  connect(m_flow_table.y, pump_m7.m_flow_in) annotation (points=[81,70; 88,70;
        88,-104; 100.7,-104], style(color=74, rgbcolor={0,0,127}));
  connect(m_flow_table.y, pump_m8.m_flow_in) annotation (points=[81,70; 88,70;
        88,-134; 100.7,-134], style(color=74, rgbcolor={0,0,127}));
end TestPressureLoss;
