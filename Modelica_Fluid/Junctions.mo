within Modelica_Fluid;
package Junctions "Junction components"
  extends Modelica_Fluid.Icons.VariantLibrary;

  model IdealTJunction
    "Splitting/joining component with static balances for an infinitesimal control volume"
    extends BaseClasses.PartialTJunction;

  equation
    connect(port_1, port_2) annotation (Line(
        points={{-100,0},{100,0}},
        color={0,127,255},
        smooth=Smooth.None));
    connect(port_1, port_3) annotation (Line(
        points={{-100,0},{0,0},{0,100}},
        color={0,127,255},
        smooth=Smooth.None));

    annotation(Documentation(info="<html>
  This model is the simplest implementation for a splitting/joining component for
  three flows. Its use is not required. It just formulates the balance
  equations in the same way that the connect symmantics would formulate them anyways.
  The main advantage of using this component is, that the user does not get
  confused when looking at the specific enthalpy at each port which might be confusing
  when not using a splitting/joining component. The reason for the confusion is that one exmanins the mixing
  enthalpy of the infinitesimal control volume introduced with the connect statement when
  looking at the specific enthalpy in the connector which
  might not be equal to the specific enthalpy at the port in the \"real world\".</html>"),
      Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));
  end IdealTJunction;

  annotation (Documentation(info="<html>
 
</html>"));
  model TJunctionVolume
    "Splitting/joining component with static balances for a dynamic control volume"
    extends BaseClasses.PartialTJunction;
    extends Volumes.BaseClasses.PartialLumpedVolume(Qs_flow = 0);

    parameter SI.Volume V "Mixing volume inside junction";

  equation
    // Only one connection allowed to a port to avoid unwanted ideal mixing
    assert(cardinality(port_1) <= 1,"
port_1 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
    assert(cardinality(port_2) <= 1,"
port_2 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
    assert(cardinality(port_3) <= 1,"
port_3 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");

    // Boundary conditions
    port_1.h_outflow = medium.h;
    port_2.h_outflow = medium.h;
    port_3.h_outflow = medium.h;

    port_1.Xi_outflow = medium.Xi;
    port_2.Xi_outflow = medium.Xi;
    port_3.Xi_outflow = medium.Xi;

    port_1.C_outflow = C;
    port_2.C_outflow = C;
    port_3.C_outflow = C;

    // Mass balances
    fluidVolume = V;
    ms_flow = port_1.m_flow + port_2.m_flow + port_3.m_flow "Mass balance";
    msXi_flow = port_1.m_flow*actualStream(port_1.Xi_outflow)
                + port_2.m_flow*actualStream(port_2.Xi_outflow)
                + port_3.m_flow*actualStream(port_3.Xi_outflow)
      "Component mass balances";

    msC_flow  = port_1.m_flow*actualStream(port_1.C_outflow)
              + port_2.m_flow*actualStream(port_2.C_outflow)
              + port_3.m_flow*actualStream(port_3.C_outflow)
      "Trace substance mass balances";

    // Momentum balance (suitable for compressible media)
    port_1.p = medium.p;
    port_2.p = medium.p;
    port_3.p = medium.p;

    // Energy balance
    Hs_flow = port_1.m_flow*actualStream(port_1.h_outflow)
              + port_2.m_flow*actualStream(port_2.h_outflow)
              + port_3.m_flow*actualStream(port_3.h_outflow);
    Ws_flow = 0;

    annotation (Documentation(info="<html>
  This model introduces a mixing volume into a junction. 
  This might be useful to examine the non-ideal mixing taking place in a real junction.</html>"),
  Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={Ellipse(
            extent={{-9,10},{11,-10}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics));
  end TJunctionVolume;

  model MultiPort
    "Multiply a port; useful if multiple connections shall be made to a port exposing a state"

    function positiveMax
      input Real x;
      output Real y;
    algorithm
      y :=max(x, 1e-10);
    end positiveMax;

    import Modelica.Constants;

    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);

    // Ports
    parameter Integer nPorts_b=1
      "Number of outlet ports (mass is distributed evenly between the outlet ports";
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
      redeclare package Medium=Medium) 
      annotation (Placement(transformation(extent={{-50,-10},{-30,10}},
            rotation=0)));
    Modelica_Fluid.Interfaces.FluidPorts_b[nPorts_b] ports_b(
      redeclare each package Medium=Medium) 
      annotation (Placement(transformation(extent={{30,40},{50,-40}},
                                  rotation=0)));

    annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-40,
              -100},{40,100}}), graphics={
          Line(
            points={{-40,0},{40,0}},
            color={0,128,255},
            thickness=1),
          Line(
            points={{-40,0},{40,26}},
            color={0,128,255},
            thickness=1),
          Line(
            points={{-40,0},{40,-26}},
            color={0,128,255},
            thickness=1),
          Text(
            extent={{-150,100},{150,60}},
            lineColor={0,0,255},
            textString="%name")}),
                            Documentation(info="<html>
<p>
This model is useful if multiple connections shall be made to a port of a volume model exposing a state,
like a pipe with ModelStructure avb. 
The mixing is shifted into the volume connected to port_a and the result is propageted back to each ports_b.
</p>
<p>
If multiple connections were directly made to the volume,
then ideal mixing would take place in the connection set, outside the volume. This is normally not intended.
</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-40,-100},{
              40,100}}),
              graphics));

  equation
    // Only one connection allowed to a port to avoid unwanted ideal mixing
    for i in 1:nPorts_b loop
      assert(cardinality(ports_b[i]) <= 1,"
each ports_b[i] of boundary shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts_b to add an additional port.
");
    end for;

    // mass and momentum balance
    0 = port_a.m_flow + sum(ports_b.m_flow);
    ports_b.p = fill(port_a.p, nPorts_b);

    // expose stream values from port_a to ports_b
    ports_b.h_outflow = fill(inStream(port_a.h_outflow), nPorts_b);
    ports_b.Xi_outflow = fill(inStream(port_a.Xi_outflow), nPorts_b);
    ports_b.C_outflow = fill(inStream(port_a.C_outflow), nPorts_b);

    // mixing at port_a
    port_a.h_outflow = sum({positiveMax(ports_b[j].m_flow)*inStream(ports_b[j].h_outflow) for j in 1:nPorts_b})
                         / sum({positiveMax(ports_b[j].m_flow) for j in 1:nPorts_b});
    for i in 1:Medium.nXi loop
      port_a.Xi_outflow[i] = sum({positiveMax(ports_b[j].m_flow)*inStream(ports_b[j].Xi_outflow[i]) for j in 1:nPorts_b})
                           / sum({positiveMax(ports_b[j].m_flow) for j in 1:nPorts_b});
    end for;
    for i in 1:Medium.nC loop
      port_a.C_outflow[i] = sum({positiveMax(ports_b[j].m_flow)*inStream(ports_b[j].C_outflow[i]) for j in 1:nPorts_b})
                           / sum({positiveMax(ports_b[j].m_flow) for j in 1:nPorts_b});
    end for;
  end MultiPort;

  package BaseClasses "Base classes for junctions"
    extends Modelica_Fluid.Icons.BaseClassLibrary;

    partial model PartialTJunction
      "Base class for a splitting/joining component with three ports"
      import Modelica_Fluid.Types;
      import Modelica_Fluid.Types.PortFlowDirection;

      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
        "Fluid medium model" 
        annotation (choicesAllMatching=true);

      Modelica_Fluid.Interfaces.FluidPort_a port_1(redeclare package Medium = 
            Medium, m_flow(min=if (portFlowDirection_1 == PortFlowDirection.Entering) then 
                    0.0 else -Modelica.Constants.inf, max=if (portFlowDirection_1
               == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
              rotation=0)));
      Modelica_Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = 
            Medium, m_flow(min=if (portFlowDirection_2 == PortFlowDirection.Entering) then 
                    0.0 else -Modelica.Constants.inf, max=if (portFlowDirection_2
               == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
        annotation (Placement(transformation(extent={{90,-10},{110,10}}, rotation=
               0)));
      Modelica_Fluid.Interfaces.FluidPort_a port_3(
        redeclare package Medium=Medium,
        m_flow(min=if (portFlowDirection_3==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
        max=if (portFlowDirection_3==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
        annotation (Placement(transformation(extent={{-10,90},{10,110}}, rotation=
               0)));

      annotation(Icon(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Rectangle(
              extent={{-100,41},{100,-47}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={192,192,192}),
            Rectangle(
              extent={{-100,37},{100,-43}},
              lineColor={0,0,0},
              fillPattern=FillPattern.HorizontalCylinder,
              fillColor={0,127,255}),
            Rectangle(
              extent={{-34,100},{34,37}},
              lineColor={0,0,0},
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={192,192,192}),
            Rectangle(
              extent={{-30,100},{30,35}},
              lineColor={0,0,0},
              fillPattern=FillPattern.VerticalCylinder,
              fillColor={0,127,255}),
            Text(
              extent={{-150,-60},{150,-100}},
              lineColor={0,0,255},
              textString="%name")}),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics));

    protected
      parameter PortFlowDirection portFlowDirection_1=PortFlowDirection.Bidirectional
        "Flow direction for port_1" 
       annotation(Dialog(tab="Advanced"));
      parameter PortFlowDirection portFlowDirection_2=PortFlowDirection.Bidirectional
        "Flow direction for port_2" 
       annotation(Dialog(tab="Advanced"));
      parameter PortFlowDirection portFlowDirection_3=PortFlowDirection.Bidirectional
        "Flow direction for port_3" 
       annotation(Dialog(tab="Advanced"));

    end PartialTJunction;
  end BaseClasses;

  package ToBeRemoved
    model MassFlowRatio "simple flow multiplier"
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);
      parameter Integer nOutlets=1
        "Number of outlet ports (mass is distributed evenly between the outlet ports";
      Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                    redeclare package Medium=Medium) 
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
              rotation=0)));
      Modelica_Fluid.Interfaces.FluidPorts_b ports_b[nOutlets] 
                                      annotation (Placement(transformation(extent=
               {{90,-40},{110,40}}, rotation=0)));

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Line(
              points={{-80,0},{80,0}},
              color={0,128,255},
              thickness=1),
            Line(
              points={{-80,0},{80,28}},
              color={0,128,255},
              thickness=1),
            Line(
              points={{-80,0},{80,-28}},
              color={0,128,255},
              thickness=1),
            Text(
              extent={{-150,100},{150,60}},
              lineColor={0,0,255},
              textString="%name")}),
                              Documentation(info="<html>
<p>
This model describes a simple flow partitioning, which is very helpful in cases where the flow is evenly distributed to several parallel flow paths which are identical in their dimensions and boundary conditions, as e.g. in heat exchangers. Only one of the parallel pipes needs to be simulated then. All flow variables in <b>port_b[i]</b> are equal to those at <b>port_a</b> divided by <b>nOutlets</b>. All effort variables are equal at all ports.
</p>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}),
                graphics));

    equation
      port_a.h_outflow  =  sum( {inStream(ports_b[i].h_outflow)     for i in 1:nOutlets})/nOutlets;
      port_a.Xi_outflow = {sum( {inStream(ports_b[i].Xi_outflow[j]) for i in 1:nOutlets})/nOutlets for j in 1:Medium.nXi};
      port_a.C_outflow  = {sum( {inStream(ports_b[i].C_outflow[j])  for i in 1:nOutlets})/nOutlets for j in 1:Medium.nXi};

      for i in 1:nOutlets loop
         ports_b[i].h_outflow  = inStream(port_a.h_outflow);
         ports_b[i].Xi_outflow = inStream(port_a.Xi_outflow);
         ports_b[i].C_outflow  = inStream(port_a.C_outflow);

         // Momentum balance
         port_a.p = ports_b[i].p;

         // Mass balance
         ports_b[i].m_flow = -port_a.m_flow/nOutlets;
      end for;
    end MassFlowRatio;

    model HeatFlowRatio "simple heat flow multiplier"
      parameter Integer nOutlets=1
        "Number of outlet ports (heat is distributed evenly between the outlet ports";
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Line(
              points={{-80,0},{80,30}},
              color={127,0,0},
              thickness=1),
            Line(
              points={{-80,0},{80,0}},
              color={127,0,0},
              thickness=1),
            Line(
              points={{-80,0},{80,-32}},
              color={127,0,0},
              thickness=1),
            Text(
              extent={{-150,100},{150,60}},
              lineColor={0,0,255},
              textString="%name")}),
                              Documentation(info="<html>
<p>
Simple model for heat flow partitioning between the two ports. The heat flow rate in port_a is divided by parameter <b>nOutlets</b> to achieve the rate at ports port_b. All temperatures are equal. The model may be used e.g. for parallel pipes in heat exchangers.
</p>
</html>"));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort_a 
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
              rotation=0)));
      Modelica_Fluid.Interfaces.HeatPorts_b heatPorts_b[nOutlets] 
                                     annotation (Placement(transformation(extent=
                {{90,-40},{110,40}}, rotation=0)));
    equation
      for i in 1:nOutlets loop
         heatPorts_b[i].Q_flow = -heatPort_a.Q_flow/nOutlets;
         heatPorts_b[i].T      =  heatPort_a.T;
      end for;
    end HeatFlowRatio;
  end ToBeRemoved;
end Junctions;
