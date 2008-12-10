within Modelica_Fluid.Examples;
package TraceSubstances "Library demonstrating the usage of trace substances"
  extends Modelica.Icons.Library;
  model RoomCO2 "Demonstrates a room volume with CO2 accumulation"
    extends Modelica.Icons.Example;
    package Medium=Modelica.Media.Air.MoistAir(extraPropertiesNames={"CO2"});
    Modelica.Blocks.Sources.Constant C(k=0.3*1.519E-3)
      "substance concentration, raising to 1000 PPM CO2" 
      annotation (Placement(transformation(extent={{-100,-28},{-80,-8}})));
    Sources.FixedBoundary boundary4(redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{80,-20},{60,0}})));
    Sensors.TraceSubstancesOnePort traceSubstanceVolume(redeclare package
        Medium = Medium) 
      annotation (Placement(transformation(extent={{0,20},{20,40}})));
    inner System system              annotation (Placement(transformation(extent={{60,60},
              {80,80}},          rotation=0)));
    Sources.PrescribedMassFlowRate_TX boundary1(
      useTraceInput=true,
      m_flow=100/1.2/3600*5,
      redeclare package Medium = Medium,
      nPorts=2) 
      annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
    Volumes.Volume volume(
      C_start={1.519E-3},
      V=100,
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      nPorts=2) annotation (Placement(transformation(extent={{-20,0},{0,20}})));
    PressureLosses.WallFrictionAndGravity pipeFriction(
      redeclare package Medium = Medium,
      length=1,
      show_Re=true,
      diameter=0.15) 
      annotation (Placement(transformation(extent={{20,-20},{40,0}})));
    Sensors.TraceSubstancesOnePort traceSubstanceSource(redeclare package
        Medium = Medium) 
      annotation (Placement(transformation(extent={{-40,20},{-20,40}})));
  equation
    connect(C.y, boundary1.C_in[1]) annotation (Line(
        points={{-79,-18},{-60,-18}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(pipeFriction.port_b, boundary4.ports[1]) annotation (Line(
        points={{40,-10},{60,-10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume.ports[2], pipeFriction.port_a) annotation (Line(
        points={{-10,-2},{-10,0},{-6,0},{-6,-10},{20,-10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(traceSubstanceVolume.port, pipeFriction.port_a) annotation (Line(
        points={{10,20},{10,-10},{20,-10}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary1.ports[1], volume.ports[1]) annotation (Line(
        points={{-40,-8},{-10,-8},{-10,2}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(boundary1.ports[2], traceSubstanceSource.port) annotation (Line(
        points={{-40,-12},{-30,-12},{-30,20}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      experiment(StopTime=3600),
      Documentation(info="<html>
This example consists of a volume with a carbon dioxide concentration that corresponds to about 1000 PPM.
There is a fresh air stream with a carbon dioxide concentration of about 300 PPM.
The fresh air stream is such that the air exchange rate is about 5 air changes per hour.
After 1 hour of ventilation, the volume's carbon dioxide concentration is close to the 
concentration of the fresh air.
</html>"));
  end RoomCO2;

  model RoomCO2WithControls "Demonstrates a room volume with CO2 controls"
    extends Modelica.Icons.Example;
    package Medium=Modelica.Media.Air.MoistAir(extraPropertiesNames={"CO2"});
    Modelica.Blocks.Sources.Constant CAtm(k=0.3*1.519E-3)
      "Atmospheric trace substance concentration, corresponding to 300 PPM CO2"
      annotation (Placement(transformation(extent={{-100,-50},{-80,-30}})));
    Sources.FixedBoundary boundary4(redeclare package Medium = Medium) 
      annotation (Placement(transformation(extent={{92,-40},{72,-20}})));
    Sensors.TraceSubstancesOnePort traceSubstanceVolume(redeclare package
        Medium = Medium) 
      annotation (Placement(transformation(extent={{0,-2},{20,18}})));
    inner System system              annotation (Placement(transformation(extent={{60,80},
              {80,100}},         rotation=0)));
    Sources.PrescribedMassFlowRate_TX freshAir(
      useTraceInput=true,
      redeclare package Medium = Medium,
      nPorts=2,
      useFlowRateInput=true) 
      annotation (Placement(transformation(extent={{-60,-42},{-40,-22}})));
    Volumes.Volume volume(
      C_start={1.519E-3},
      V=100,
      redeclare package Medium = Medium,
      initType=Modelica_Fluid.Types.Init.InitialValues,
      nPorts=3) annotation (Placement(transformation(extent={{-20,-22},{0,-2}})));
    PressureLosses.WallFrictionAndGravity pipeFriction(
      redeclare package Medium = Medium,
      length=1,
      show_Re=true,
      diameter=0.15) 
      annotation (Placement(transformation(extent={{32,-40},{52,-20}})));
    Sensors.TraceSubstancesOnePort traceSubstanceSource(redeclare package
        Medium = Medium) 
      annotation (Placement(transformation(extent={{-40,-2},{-20,18}})));
    Sources.PrescribedMassFlowRate_TX peopleSource(
      m_flow=100/1.2/3600*5,
      redeclare package Medium = Medium,
      nPorts=1,
      useFlowRateInput=true,
      useTraceInput=false,
      C={100}) "CO2 emitted by room occupants." 
      annotation (Placement(transformation(extent={{-38,-98},{-18,-78}})));
    Modelica.Blocks.Sources.CombiTimeTable NumberOfPeople(table=[0,0; 9*3600,0;
          9*3600,10; 11*3600,10; 11*3600,2; 13*3600,2; 13*3600,15; 15*3600,15;
          15*3600,5; 18*3600,5; 18*3600,0; 24*3600,0])
      "Time table for number of people in the room" 
      annotation (Placement(transformation(extent={{-100,-90},{-80,-70}})));
    Modelica.Blocks.Math.Gain gain(k=8.18E-6/100)
      "CO2 mass flow rate, released per 100 person (there is another 100 factor in peopleSource)"
      annotation (Placement(transformation(extent={{-68,-90},{-48,-70}})));
    Modelica.Blocks.Math.Gain gain1(k=-100*1.2/3600*5)
      "Nominal fresh air flow rate (for u=1)" 
      annotation (Placement(transformation(extent={{0,40},{20,60}})));
    Modelica.Blocks.Math.Gain gainSensor(k=1/1.519E-3)
      "Gain to normalize CO2 measurement signal. y=1 corresponds to 1000 PPM" 
      annotation (Placement(transformation(extent={{40,-2},{60,18}})));
    Modelica.Blocks.Sources.Constant CO2Set(k=1) "Normalized CO2 set point" 
      annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
    Modelica.Blocks.Continuous.LimPID PID(
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      yMax=0,
      yMin=-1,
      Ti=10,
      k=10)   annotation (Placement(transformation(extent={{-40,40},{-20,60}})));
  equation
    connect(CAtm.y, freshAir.C_in[1]) 
                                    annotation (Line(
        points={{-79,-40},{-60,-40}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(pipeFriction.port_b, boundary4.ports[1]) annotation (Line(
        points={{52,-30},{72,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(volume.ports[2], pipeFriction.port_a) annotation (Line(
        points={{-10,-22},{-10,-30},{32,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(traceSubstanceVolume.port, pipeFriction.port_a) annotation (Line(
        points={{10,-2},{10,-30},{32,-30}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(freshAir.ports[1], volume.ports[1])  annotation (Line(
        points={{-40,-30},{-10,-30},{-10,-19.3333}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(freshAir.ports[2], traceSubstanceSource.port)  annotation (Line(
        points={{-40,-34},{-30,-34},{-30,-2}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics),
      experiment(StopTime=86400, tolerance=1e-006),
      Documentation(info="<html>
This example illustrates a room volume with a CO2 source and a fresh air supply with feedback
control.
The CO2 emission rate is proportional to the room occupancy, which is defined by a schedule.
The fresh air flow rate is controlled such that the room CO2
concentration does not exceed <tt>1000 PPM (=1.519E-3 kg/kg)</tt>.
The fresh air has a CO2 concentration of <tt>300 PPM</tt> which corresponds to a typical
CO2 concentration in the outside air. 
</p>
<p>
The CO2 emission from the occupants is implemented as a mass flow source.
Depending on the activity and size, a person emits about <tt>8.18E-6 kg/s</tt> CO2. In the model, 
this value is multiplied by the number of occupants. 
Since the mass flow rate associate with the CO2 source model contributes to the volume's energy balance,
this mass flow rate should be kept small. Thus, in the source model, we set the
CO2 concentration to <tt>C={100} kg/kg</tt>, and scaled the mass flow rate using
<pre>
  m_flow = 1/100 * nPeo * 8.18E-6 kg/(s*person)
</pre>
where <tt>nPeo</tt> is the number of people in the room.
This results in a mass flow rate that is about 5 orders of magnitudes smaller than the supply air flow rate,
and hence its contribution to the volume's energy balance is negligible.
</html>"));
    connect(NumberOfPeople.y[1], gain.u) annotation (Line(
        points={{-79,-80},{-70,-80}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(gain.y, peopleSource.m_flow_in) annotation (Line(
        points={{-47,-80},{-38,-80}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(peopleSource.ports[1], volume.ports[3]) annotation (Line(
        points={{-18,-88},{-8,-88},{-8,-24.6667},{-10,-24.6667}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(traceSubstanceVolume.C, gainSensor.u) 
                                             annotation (Line(
        points={{21,8},{38,8}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(CO2Set.y, PID.u_s) annotation (Line(
        points={{-59,50},{-42,50}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(gainSensor.y, PID.u_m) 
                              annotation (Line(
        points={{61,8},{70,8},{70,30},{-30,30},{-30,38}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(PID.y, gain1.u) annotation (Line(
        points={{-19,50},{-2,50}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(gain1.y, freshAir.m_flow_in)  annotation (Line(
        points={{21,50},{30,50},{30,70},{-88,70},{-88,-24},{-60,-24}},
        color={0,0,127},
        smooth=Smooth.None));
  end RoomCO2WithControls;
end TraceSubstances;
