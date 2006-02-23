package ControlValves 
    model ValveIncompressible "Valve for (almost) incompressible fluids" 
    extends BaseClasses.ControlValves.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    annotation (
    Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluids. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>This model can be used with any low compressibility fluids, such as liquids or gases at very low pressure drops.
</html>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
    initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
    
    equation 
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
               smooth(0,if dp>=0 then sqrtR(dp) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dp);
    end if;
    end ValveIncompressible;
  
  model ValveVaporizing 
    "Valve for possibly vaporizing (almost) incompressible fluids, accounts for choked flow conditions" 
    import Modelica_Fluid.Types.CvTypes;
    extends BaseClasses.ControlValves.PartialValve(
      redeclare replaceable package Medium = 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    parameter Real Fl_nom=0.9 "Liquid pressure recovery factor";
    replaceable function FlCharacteristic = 
        Modelica_Fluid.BaseClasses.ControlValves.ValveCharacteristics.one 
      extends 
      Modelica_Fluid.BaseClasses.ControlValves.ValveCharacteristics.baseFun 
      "Pressure recovery characteristic";
    Real Ff "Ff coefficient (see IEC/ISA standard)";
    Real Fl "Pressure recovery coefficient Fl (see IEC/ISA standard)";
    Medium.AbsolutePressure pv "Saturation pressure";
    SI.Pressure dpEff "Effective pressure drop";
    Medium.AbsolutePressure pin "Inlet pressure";
    Medium.AbsolutePressure pout "Outlet pressure";
    annotation (
      Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
      Diagram,
      Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid, with possible choked flow conditions. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).<p>
The model operating range includes choked flow operation, which takes place for low outlet pressures due to flashing in the vena contracta; otherwise, non-choking conditions are assumed.
<p>This model can be used with two-phase medium models, to describe the liquid and (possible) two-phase conditions.
<p>The default liquid pressure recovery coefficient <tt>Fl</tt> is constant and given by the parameter <tt>Fl_nom</tt>. The relative change (per unit) of the recovery coefficient can be specified as a given function of the valve opening by replacing the <tt>FlCharacteristic</tt> function.
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
  equation 
    pin = port_a.p;
    pout = port_b.p;
    pv = Medium.saturationPressure(T);
    Ff = 0.96 - 0.28*sqrt(pv/Medium.fluidConstants[1].criticalPressure);
    Fl = Fl_nom*FlCharacteristic(stemPosition);
    dpEff = if pout < (1 - Fl^2)*pin + Ff*Fl^2*pv then 
              Fl^2*(pin - Ff*pv) else dp 
      "Effective pressure drop, accounting for possible choked conditions";
    if CheckValve then
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
           (if dpEff>=0 then sqrtR(dpEff) else 0);
     else
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dpEff);
    end if;
  end ValveVaporizing;
  
  model ValveCompressible 
    "Valve for compressible fluids, accounts for choked flow conditions" 
    extends BaseClasses.ControlValves.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    parameter Real Fxt_full=0.5 "Fk*xt critical ratio at full opening";
    replaceable function xtCharacteristic = 
        Modelica_Fluid.BaseClasses.ControlValves.ValveCharacteristics.one 
      extends 
      Modelica_Fluid.BaseClasses.ControlValves.ValveCharacteristics.baseFun 
      "Critical ratio characteristic";
    Real Fxt;
    Real x "Pressure drop ratio";
    Real xs "Saturated pressure drop ratio";
    Real Y "Compressibility factor";
    Medium.AbsolutePressure p "Inlet pressure";
  protected 
    parameter Real Fxt_nom(fixed=false) "Nominal Fxt";
    parameter Real x_nom(fixed=false) "Nominal pressure drop ratio";
    parameter Real xs_nom(fixed=false) "Nominal saturated pressure drop ratio";
    parameter Real Y_nom(fixed=false) "Nominal compressibility factor";
    
    annotation (
    Icon(Text(extent=[-100, -40; 100, -80], string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, compressible fluid. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>The product Fk*xt is given by the parameter <tt>Fxt_full</tt>, and is assumed constant by default. The relative change (per unit) of the xt coefficient with the valve opening can be specified by replacing the <tt>xtCharacteristic</tt> function.
</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      // Determination of Av by the nominal operating point conditions
      Fxt_nom = Fxt_full*xtCharacteristic(stemPosition_nom);
      x_nom = dp_nom/p_nom;
      xs_nom = smooth(0, if x_nom > Fxt_nom then Fxt_nom else x_nom);
      Y_nom = 1 - abs(xs_nom)/(3*Fxt_nom);
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*Y_nom*sqrt(d_nom)*sqrtR(p_nom*xs_nom);
    else
      // Dummy values
      Fxt_nom = 0;
      x_nom = 0;
      xs_nom = 0;
      Y_nom = 0;
    end if;
    
  equation 
    p = noEvent(if dp>=0 then port_a.p else port_b.p);
    Fxt = Fxt_full*xtCharacteristic(stemPosition);
    x = dp/p;
    xs = smooth(0, if x < -Fxt then -Fxt else if x > Fxt then Fxt else x);
    Y = 1 - abs(xs)/(3*Fxt);
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*
        smooth(0,if xs>=0 then sqrtR(p*xs) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*sqrtR(p*xs);
    end if;
  end ValveCompressible;
  
  model ValveLinear "Valve for water/steam flows with linear pressure drop" 
    extends BaseClasses.Common.PartialTwoPortTransport;
    parameter Types.HydraulicConductance Kv 
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20,70; 20,110],   rotation=-90);
  equation 
    m_flow = Kv*opening*dp;
    
  annotation (
    Icon(
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
           Text(extent=[-143,-66; 148,-106],  string="%name")),
    Diagram,
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
</html>"),
      Coordsys(grid=[1,1], scale=0));
  end ValveLinear;
  
  model ValveDiscrete "Valve for water/steam flows with linear pressure drop" 
    extends Modelica_Fluid.BaseClasses.Common.PartialTwoPortTransport;
    parameter Modelica_Fluid.Types.HydraulicConductance Kv 
      "Hydraulic conductance for open valve (m_flow = Kv*dp)";
    parameter Real Kv_small_rel = 0 
      "Relative hydraulic conductance for closed valve (m_flow = Kv_small_rel*Kv*dp)";
    Modelica.Blocks.Interfaces.BooleanInput open 
    annotation (extent=[-20,60; 20,100],   rotation=-90);
  equation 
    m_flow = if open then Kv*dp else Kv_small_rel*Kv*dp;
    
  annotation (
    Icon(
        Line(points=[0,50; 0,0],   style(
            color=0,
            rgbcolor={0,0,0},
            fillPattern=1)),
        Rectangle(extent=[-20,60; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
           Text(extent=[-145,-58; 146,-98],   string="%name"),
        Polygon(points=[-100,50; 100,-50; 100,50; 0,0; -100,-50; -100,50], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=DynamicSelect(7, if open > 0.5 then 2 else 7)))),
    Diagram,
    Documentation(info="<HTML>
<
<p>
This very simple model provides a pressure drop which is proportional to the flowrate if the Boolean open signal is <b>true</b>. Otherwise, the
mass flow rate is zero. If Kv_small_rel > 0, a small leakage
mass flow rate occurs when open = <b>false</b>. This might be
useful in certain situations when the model is not
mathematically well-defined due to a closed valve.
</p>
<p>
In a diagram animation, the valve is shown in \"green\", when
it is open.
</p>
</HTML>",
      revisions="<html>
<ul>
<li><i>Nov 2005</i>
    by Katja Poschlad (based on ValveLinear).</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
  end ValveDiscrete;
  
end ControlValves;
