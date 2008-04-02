within FluidSandbox.Examples.StandardConnectionSemantics;
model NonlinearTest 
  extends Icons.Example;
  import Modelica.SIunits.Conversions;
  
  annotation (
    experiment(StopTime=4, Tolerance=1e-006),
    experimentSetupOutput,
    Documentation(info="<html>
<body>
<h2>Exercise 4 of the Modelica Fluid Tutorial</h2>
<h3>Non-linear systems of equations</h3>
<p>In this model a few pressure drop models are combined with one junction point and connected directly.
This is possible, but can have drawbacks.</p>
 
<h4>Exercises</h4>
Run the model as it is given in Exercise 4. Repeat with other liquids, e.g. water with states p and h or p and T, or a mixture of ideal gases. What are the sizes of the non-linear equation systems? Why do they differ? Look out both for the equation system in the dynamic part and in the initialization part.<p>
<p>With the given parameters the non-linear equation system in  the model is difficult to solve. This will be demonstrated with the following exercises</p>
<h4>Exercises</h4>
Run the model as it is given in Exercise 4 in the package ThermodynamicTutorial, for 4 seconds
(settings are stored). For all runs and settings, observe the sizes of the non-linear equation systems that occur.
<ol>
<li>Run the model with a compressible medium from WaterIF97, using p and h as states. The model will fail. Re-run with initialization set to InitialStates. Now it works, with a smaller system of equations.</li>
<li>Duplicate the mixing volume model and insert it after the pipe component. Re-run with different settings for the initialization and observe the sizes of the initialization system.</li>
<li>Choose a the medium \"Simple natural gas mixture with 6 components\" as medium and run the exercise with and without an additional mixing volume, and with steady state and initial value initial settings.</li>
</ol>
<p><b>Revision history:</b></p>
<ul>
<li><i>March 3 2005</i>
    by <a href=\"mailto:hubertus@modelon.se\">Hubertus Tummescheit</a>:<br>
<li><i>August 22 2006</i>
    updated to latest MiniFluid version and extended, <a href=\"mailto:hubertus@modelon.se\">Hubertus Tummescheit</a>:<br>
</body>
</html>"),
    Diagram);
  Sources.PrescribedBoundary_pTX_B Sink(
    p=500000, 
    T=Conversions.from_degC(20), 
    redeclare package FluidInterface = FluidInterface, 
    redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
    annotation (extent=[100,0; 80,20], rotation=0);
  Sources.PrescribedMassFlowRate_TX_B MassFlowSource(
    useFlowRateInput=true, 
    T=Conversions.from_degC(20), 
    redeclare package FluidInterface = FluidInterface, 
    redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
    annotation (extent=[-80,0; -60,20], rotation=
            0);
  Modelica.Blocks.Sources.Ramp flowsource(
    duration=2,
    startTime=1,
    offset=1,
    height=1)   annotation (extent=[-100,46; -80,
            66],
          rotation=0);
  PressureLosses.WallFriction pipe(
    roughness=0.0002,
    length=1,
    diameter=0.1,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false, 
    redeclare package FluidInterface = FluidInterface, 
    redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
                   annotation (extent=[10,0; 30,20],
          rotation=0);
  PressureLosses.WallFrictionAA pipe1(
    roughness=0.0002,
    length=1,
    diameter=0.1,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false, 
    redeclare package FluidInterface = FluidInterface, 
    redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
                   annotation (extent=[50,0; 70,20],
          rotation=0);
   PressureLosses.WallFriction pipe2(
    roughness=0.0002,
    length=1,
    diameter=0.1,
    provide_p_a=false,
    provide_p_b=false,
    provide_T_a=false,
    provide_T_b=false,
    provide_m_flow_ab=false, 
    redeclare package FluidInterface = FluidInterface, 
    redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
                   annotation (extent=[40,-40; 60,-20],
                  rotation=0);
  Sources.PrescribedBoundary_pTX_A Sink1(
    p=510000, 
    T=Conversions.from_degC(20), 
    redeclare package FluidInterface = FluidInterface, 
    redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
    annotation (extent=[100,-40; 80,-20],
          rotation=0);
  Volumes.Volume mixingVolume(
    V=1e-5,
    provide_p=false,
    provide_T=false, 
    initType=Modelica_Fluid.Types.Init.InitialValues, 
    redeclare package FluidInterface = FluidInterface, 
    p_start=501000, 
    redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
                    annotation (extent=[-50,0; -30,20],
                 rotation=0);
  Junctions.IdealJunctionAAB idealJunctionAAB(redeclare package FluidInterface 
      = FluidInterface, redeclare package Medium = 
        Modelica.Media.IdealGases.MixtureGases.FlueGasSixComponents) 
    annotation (extent=[-20,0; 0,20], rotation=180);
equation 
  
  connect(MassFlowSource.port, mixingVolume.port_a[1]) annotation (points=[-60,10; 
        -50,10],     style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(idealJunctionAAB.port_1, pipe.port_a) annotation (points=[0,10; 5,10; 
        5,10; 10,10], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(pipe1.port_a, pipe.port_b) annotation (points=[50,10; 30,10], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(idealJunctionAAB.port_3, pipe2.port_a) annotation (points=[-10,0; -10,
        -30; 40,-30], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(pipe2.port_b, Sink1.port) annotation (points=[60,-30; 80,-30], style(
      color=69,
      rgbcolor={0,127,255},
      smooth=0));
  connect(flowsource.y, MassFlowSource.m_flow_in) annotation (points=[-79,56;
        -70,56; -70,30; -90,30; -90,16; -79.3,16], style(
      color=74,
      rgbcolor={0,0,127},
      smooth=0));
  connect(pipe1.port_b, Sink.port) annotation (points=[70,10; 80,10], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
  connect(mixingVolume.port_b[1], idealJunctionAAB.port_2) annotation (points=[
        -30,10; -25,10; -25,10; -20,10], style(
      color=69, 
      rgbcolor={0,127,255}, 
      smooth=0));
end NonlinearTest;
