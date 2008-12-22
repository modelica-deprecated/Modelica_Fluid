within Modelica_Fluid;
package Interfaces
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow"

  annotation (Documentation(info="<html>
 
</html>", revisions="<html>
<ul>
<li><i>June 9th, 2008</i>
       by Michael Sielemann: Introduced stream keyword after decision at 57th Design Meeting (Lund).</li>
<li><i>May 30, 2007</i>
       by Christoph Richter: moved everything back to its original position in Modelica_Fluid.</li>
<li><i>Apr. 20, 2007</i>
       by Christoph Richter: moved parts of the original package from Modelica_Fluid
       to the development branch of Modelica 2.2.2.</li>
<li><i>Nov. 2, 2005</i>
       by Francesco Casella: restructured after 45th Design Meeting.</li>
<li><i>Nov. 20-21, 2002</i>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.
<li><i>Nov. 11, 2002</i>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><i>Nov. 6, 2002</i>
       by Hilding Elmqvist: first version.</li>
<li><i>Aug. 11, 2002</i>
       by Martin Otter: Improved according to discussion with Hilding
       Elmqvist and Hubertus Tummescheit.<br>
       The PortVicinity model is manually
       expanded in the base models.<br>
       The Volume used for components is renamed
       PartialComponentVolume.<br>
       A new volume model \"Fluid.Components.PortVolume\"
       introduced that has the medium properties of the port to which it is
       connected.<br>
       Fluid.Interfaces.PartialTwoPortTransport is a component
       for elementary two port transport elements, whereas PartialTwoPort
       is a component for a container component.</li>
</li>
</ul>
</html>"));

  extends Modelica.Icons.Library;

  connector FluidPort
    "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)"

    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
      "Medium model" annotation (choicesAllMatching=true);

    flow Medium.MassFlowRate m_flow
      "Mass flow rate from the connection point into the component";
    Medium.AbsolutePressure p "Pressure in the connection point";
    stream Medium.SpecificEnthalpy h_outflow
      "Specific enthalpy close to the connection point if m_flow < 0";
    stream Medium.MassFraction Xi_outflow[Medium.nXi]
      "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";
    stream Medium.ExtraProperty C_outflow[Medium.nC]
      "Properties c_i/m close to the connection point if m_flow < 0";
  end FluidPort;

  connector FluidPort_a "Generic fluid connector at design inlet"
    extends FluidPort;
    annotation (defaultComponentName="port_a",
                Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid), Text(extent={{-150,110},{150,50}},
              textString="%name")}),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,127,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid), Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPort_a;

  connector FluidPort_b "Generic fluid connector at design outlet"
    extends FluidPort;
    annotation (defaultComponentName="port_b",
                Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics={
          Ellipse(
            extent={{-40,40},{40,-40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(extent={{-150,110},{150,50}}, textString="%name")}),
         Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics={
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,127,255},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPort_b;

  connector FluidPorts_a
    "Fluid connector with filled, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
    extends FluidPort;
    annotation (defaultComponentName="ports_a",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString="%name"),
          Rectangle(
            extent={{25,-100},{-25,100}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,90},{25,40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,25},{25,-25}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,-40},{25,-90}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{50,-200},{-50,200}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,180},{50,80}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,-80},{50,-180}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPorts_a;

  connector FluidPorts_b
    "Fluid connector with outlined, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
    extends FluidPort;
    annotation (defaultComponentName="ports_b",
                Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Text(extent={{-75,130},{75,100}}, textString="%name"),
          Rectangle(
            extent={{-25,100},{25,-100}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,90},{25,40}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,25},{25,-25}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-25,-40},{25,-90}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,-50},{15,-80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,15},{15,-15}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-15,50},{15,80}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}),
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-50,-200},{50,200}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-50,200},{50,-200}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,180},{50,80}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,50},{50,-50}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-50,-80},{50,-180}},
            lineColor={0,0,0},
            fillColor={0,127,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,30},{30,-30}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,100},{30,160}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-30,-100},{30,-160}},
            lineColor={0,127,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end FluidPorts_b;

  connector HeatPorts_a
    "HeatPort connector with filled, large icon to be used for vectors of HeatPorts (vector dimensions must be added after dragging)"
    extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
    annotation (defaultComponentName="heatPorts_a",
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-200,-50},{200,50}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-201,50},{200,-50}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-171,45},{-83,-45}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-45,45},{43,-45}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{82,45},{170,-45}},
            lineColor={127,0,0},
            fillColor={127,0,0},
            fillPattern=FillPattern.Solid)}));
  end HeatPorts_a;

  connector HeatPorts_b
    "HeatPort connector with filled, large icon to be used for vectors of HeatPorts (vector dimensions must be added after dragging)"
    extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
    annotation (defaultComponentName="heatPorts_b",
         Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-200,-50},{200,50}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-200,50},{200,-51}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-170,44},{-82,-46}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-44,46},{44,-44}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{82,45},{170,-45}},
            lineColor={127,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid)}));
  end HeatPorts_b;

  partial model PartialTwoPort "Partial component with two ports"
    import Modelica.Constants;
    outer Modelica_Fluid.System system "System wide properties";

    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    parameter Boolean allowFlowReversal = system.allowFlowReversal
      "= true to allow flow reversal, false restricts to design direction (port_a -> port_b)"
      annotation(Dialog(tab="Assumptions"), Evaluate=true);

    Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                  redeclare package Medium = Medium,
                       m_flow(min=if allowFlowReversal then -Constants.inf else 0))
      "Fluid connector a (positive design flow direction is from port_a to port_b)"
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}},
              rotation=0)));
    Modelica_Fluid.Interfaces.FluidPort_b port_b(
                                  redeclare package Medium = Medium,
                       m_flow(max=if allowFlowReversal then +Constants.inf else 0))
      "Fluid connector b (positive design flow direction is from port_a to port_b)"
      annotation (Placement(transformation(extent={{110,-10},{90,10}}, rotation=
               0), iconTransformation(extent={{110,-10},{90,10}})));
    // Model structure, e.g. used for visualization
  protected
    parameter Boolean port_a_exposesState = false
      "= true if port_a exposes the state of a fluid volume";
    parameter Boolean port_b_exposesState = false
      "= true if port_b.p exposes the state of a fluid volume";

    annotation (
      Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics),
      Documentation(info="<html>
<p>
This partial model defines an interface for components with two ports. 
The components may transport fluid and have internal storage. 
The treatment of design flow direction and flow reversal are predefined. The variable
</p>
<ul>
<li><tt>m_flow</tt> defines the mass flow rate in design direction, which is port_a.m_flow.</li> 
</ul>
<p>
A derived component providing direct access to internal storage of mass or energy through port_a or port_b 
should redefine the protected parameters port_a_exposesState and port_b_exposesState appropriately. 
This will be visualized at the port icons, in order to improve the understanding of fluid model diagrams.
</p>
</html>"),
      Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
          Polygon(
            points={{20,-70},{60,-85},{20,-100},{20,-70}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{20,-75},{50,-85},{20,-95},{20,-75}},
            lineColor={255,255,255},
            smooth=Smooth.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            visible=allowFlowReversal),
          Line(
            points={{55,-85},{-60,-85}},
            color={0,128,255},
            smooth=Smooth.None),
          Text(
            extent={{-149,-114},{151,-154}},
            lineColor={0,0,255},
            fillPattern=FillPattern.HorizontalCylinder,
            fillColor={0,127,255},
            textString="%name"),
          Ellipse(
            extent={{-110,26},{-90,-24}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            visible=port_a_exposesState),
          Ellipse(
            extent={{90,25},{110,-25}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid,
            visible=port_b_exposesState)}));
  end PartialTwoPort;

      partial model PartialPressureLoss
    "Common interface for pressure loss models"

        // Medium
        replaceable package Medium = 
          Modelica.Media.Interfaces.PartialMedium "fluid medium" 
            annotation(Dialog(tab="Internal Interface", enable=false));

        // Discretization
        parameter Integer n=1 "number of flow segments" 
          annotation(Dialog(tab="Internal Interface", enable=false));

        input Medium.ThermodynamicState[n+1] state "states along design flow" 
           annotation(Dialog(tab="Internal Interface", enable=false,group="Inputs"));
        output Medium.MassFlowRate[n] m_flow
      "mass flow rates along design flow" 
           annotation(Dialog(tab="Internal Interface", enable=false,group="Outputs"));

        // Variables
        Medium.AbsolutePressure[n+1] p "pressures of states";
        Modelica.SIunits.Pressure[n] dp "pressure loss between states";

      equation
        p = Medium.pressure(state);
        dp = p[1:n] - p[2:n+1];

        annotation (
           Documentation(info="<html>
<p>
The mass flow rates <tt>m_flow[n]</tt> between n+1 flow segments 
are obtained as function of the thermodynamic <tt>state</tt> of the flow segments for a given fluid <tt>Medium</tt>.
</p>
<p>
An extending model implementing this interface needs to define the relation between the predefined pressure drops <tt>dp[n]</tt> 
and <tt>m_flow[n]</tt>.
</p>
</html>"), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
              {100,100}}), graphics={Line(
            points={{-80,-50},{-80,50},{80,-50},{80,50}},
            color={0,0,255},
            smooth=Smooth.None,
            thickness=1)}));

      end PartialPressureLoss;

  partial model PartialHeatTransfer "Common interface for heat transfer models"

    // Parameters
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium 
      annotation(Dialog(tab="Internal Interface", enable=false));

    parameter Integer n=1 "Number of heat transfer segments" 
      annotation(Dialog(tab="Internal Interface", enable=false), Evaluate=true);

    // Inputs provided to heat transfer model
    input Medium.ThermodynamicState[n] state 
      annotation(Dialog(tab="Internal Interface", enable=false));

    // Output defined by heat transfer model
    output SI.HeatFlowRate[n] Q_flow "Heat flow rates per tube";

    // Heat ports
    Modelica_Fluid.Interfaces.HeatPorts_a[n] heatPorts
      "Heat port to component boundary" 
      annotation (Placement(transformation(extent={{-10,60},{10,80}},
              rotation=0), iconTransformation(extent={{-20,60},{20,80}})));

    // Variables
    SI.Temperature[n] T;

  equation
    T = Medium.temperature(state);
    heatPorts.Q_flow = Q_flow;

    annotation (Documentation(info="<html>
<p>
The heat flow rates <tt>Q_flow[n]</tt> through the boundaries of n flow segments 
are obtained as function of the thermodynamic <tt>state</tt> of the flow segments for a given fluid <tt>Medium</tt>
and the boundary temperatures <tt>heatPorts[n].T</tt>.
</p>
<p>
An extending model implementing this interface needs to define the relation between the predefined fluid temperatures <tt>T[n]</tt>,
the boundary temperatures <tt>heatPorts[n].T</tt>, and the heat flow rates <tt>Q_flow[n]</tt>.
</p>
</html>"),
        Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-60,64},{60,-56}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={232,0,0}), Text(
            extent={{-38,26},{40,-14}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={232,0,0},
            textString="%name")}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}),
                      graphics));
  end PartialHeatTransfer;
end Interfaces;
