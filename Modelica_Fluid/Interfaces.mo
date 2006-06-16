package Interfaces 
  "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow" 
  
  annotation (Documentation(info="<html>
</html>", revisions="<html>
<ul>
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
  import SI = Modelica.SIunits;
  
  
  
  
  package ControlVolumes 
  partial model PartialTwoPortTransport 
      "Partial element transporting fluid between two ports without storing mass or energy" 
      import SI = Modelica.SIunits;
      import Modelica.Constants;
    replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
                                                                         annotation (
        choicesAllMatching =                                                                            true);
      
   parameter Types.FlowDirectionWithGlobalDefault.Temp flowDirection=
                     Modelica_Fluid.Types.FlowDirection.Unidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
      
    Modelica_Fluid.Interfaces.Ports.FluidPort_a port_a(
                                  redeclare package Medium = Medium,
                       m_flow(start=0,min=if allowFlowReversal then -Constants.inf else 0)) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
      annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.Ports.FluidPort_b port_b(
                                  redeclare package Medium = Medium,
                       m_flow(start=0,max=if allowFlowReversal then +Constants.inf else 0)) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
      annotation (extent=[110,-10; 90,10]);
    Medium.BaseProperties medium_a "Medium properties in port_a";
    Medium.BaseProperties medium_b "Medium properties in port_b";
    Medium.MassFlowRate m_flow(start=0) 
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
    SI.VolumeFlowRate V_flow_a = port_a.m_flow/medium_a.d 
        "Volume flow rate near port_a";
    SI.Pressure dp(start=0) "Pressure difference between port_a and port_b";
      
    annotation (
      Coordsys(grid=[1, 1], component=[20, 20]),
      Diagram,
      Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"),
      Icon);
    protected 
      parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
  equation 
    // Properties in the ports
    port_a.p   = medium_a.p;
    port_a.h   = medium_a.h;
    port_a.Xi = medium_a.Xi;
    port_b.p   = medium_b.p;
    port_b.h   = medium_b.h;
    port_b.Xi = medium_b.Xi;
      
    /* Handle reverse and zero flow */
    port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,  port_b.h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
      
    /* Energy, mass and substance mass balance */
    port_a.H_flow + port_b.H_flow = 0;
    port_a.m_flow + port_b.m_flow = 0;
    port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
      
    // Design direction of mass flow rate
    m_flow = port_a.m_flow;
      
    // Pressure difference between ports
    dp = port_a.p - port_b.p;
      
  end PartialTwoPortTransport;
    
  partial model PartialDistributedFlow 
      import Modelica_Fluid.Types;
      import Modelica.Constants.*;
    replaceable package Medium = PackageMedium 
      extends Modelica.Media.Interfaces.PartialMedium "Fluid medium model" 
        annotation (choicesAllMatching=true);
      
  //Discretization
    parameter Integer n(min=1)=1 "Number of pipe segments";
    final parameter Integer np=if singleState_hydraulic then 2 else n + 1 
        "Number of momentum balances"                                                     annotation(Dialog(tab="Advanced"),Evaluate=true);
    final parameter Integer nl=integer(n/2)+1 
        "Number of control volume that contains single state"                 annotation(Evaluate=true);
      
  //Advanced model options
    parameter Boolean allowFlowReversal=true 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
                                                                       annotation(Dialog(tab="Advanced", group="Mass and energy balances", enable=not static));
    parameter Boolean static=false 
        "= true, static balances, no mass or energy is stored" 
                                  annotation(Dialog(tab="Advanced", group="Mass and energy balances"),Evaluate=true);
    parameter Boolean singleState_hydraulic=false 
        " = true, lumped pressure drop, reduces number of pressure states to one"
                                                                                annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
    parameter Boolean singleState_thermal=false 
        " = true, number of temperature or enthalpy states is reduced to one" 
                                                                             annotation(Evaluate=true,Dialog(tab="Advanced", group="Mass and energy balances", enable=not static));
      
  //Initialization
      parameter Types.Init.Temp initType=Types.
          Init.NoInit "Initialization option" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Boolean use_T_start=true 
        "Use T_start if true, otherwise h_start" 
      annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_a_start=Medium.p_default 
        "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start=Medium.p_default 
        "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));
      final parameter Medium.AbsolutePressure[n] p_start=if n > 1 then linspace(
          p_a_start - (p_a_start - p_b_start)/(2*n),
          p_b_start + (p_a_start - p_b_start)/(2*n),
          n) else {(p_a_start + p_b_start)/2} "Start value of pressure";
      parameter Medium.Temperature T_start=if use_T_start then Medium.T_default else 
                Medium.temperature_phX(
          (p_a_start + p_b_start)/2,
          h_start,
          X_start) "Start value of temperature" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=if use_T_start then 
          Medium.specificEnthalpy_pTX(
          (p_a_start + p_b_start)/2,
          T_start,
          X_start) else Medium.h_default "Start value of specific enthalpy" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
        "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
      parameter Medium.MassFlowRate mflow_start 
        "Start value for mass flow rate"                                       annotation(Evaluate=true, Dialog(tab = "Initialization"));
      final parameter SI.Pressure dp_start=p_a_start - p_b_start;
      
  //Geometry parameters
    parameter Modelica_Fluid.Types.CrossSectionTypes.Temp crossSectionType=
        Modelica_Fluid.Types.CrossSectionTypes.Circular 
        "Type of cross section of pipe" 
      annotation (Dialog(tab="General", group="Geometry"));
    parameter SI.Length length "Length"   annotation(Dialog(tab="General", group="Geometry"));
    parameter SI.Diameter diameter "Diameter of circular pipe"      annotation(Dialog(group="Geometry", enable=crossSectionType==1));
    parameter SI.Length height "Height of rectangular pipe"           annotation(Dialog(group="Geometry", enable=crossSectionType==2));
    parameter SI.Length width "Width of rectangular pipe"            annotation(Dialog(group="Geometry", enable=crossSectionType==2));
    parameter SI.Length perimeter=if crossSectionType == 1 then Modelica.Constants.pi*diameter else if crossSectionType == 2 then 2*height + 2*
        width else 1 "Inner perimeter"  annotation(Dialog(tab="General", group="Geometry", enable=crossSectionType==3));
    inner parameter SI.Area area=if crossSectionType == 1 then Modelica.Constants.pi*diameter*diameter/4 else if crossSectionType
         == 2 then height*width else 1 "Inner cross section area" 
                                            annotation(Dialog(tab="General", group="Geometry", enable=crossSectionType==3));
    SI.Volume[n] Vi "Discretized volume, determine in inheriting class ";
      
  //Total quantities
    SI.Energy[n] U "Internal energy of fluid";
    SI.Mass[n] m "Mass of fluid";
    SI.Mass[n,Medium.nXi] mXi "Masses of independent components in the fluid";
      
  //Flow quantities
    inner Medium.MassFlowRate[n + 1] m_flow(each min=if allowFlowReversal then -inf else 
                0, each start=mflow_start, each fixed=false) 
        "Mass flow rates of fluid across segment boundaries";
    SI.Velocity[n+1] v "velocity at volume boundaries (for display purposes)";
    Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow 
        "Independent mass flow rates across segment boundaries";
    Medium.EnthalpyFlowRate[n + 1] H_flow 
        "Enthalpy flow rates of fluid across segment boundaries";
   parameter Boolean use_d_nominal=false 
        "= true, if d_nominal is used, otherwise computed from medium"                              annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
   parameter SI.Density d_nominal=0.01 
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
                                                               annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_nominal));
      
  //Source terms, have to be set in inheriting class
    Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
    Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
        "Independent mass flow rates, source or sink";
    SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
    SI.Pressure[np] dp(start=dp0) "pressure difference across staggered grid";
      
  //Fluid ports
    Modelica_Fluid.Interfaces.Ports.FluidPort_a port_a(
                                                 redeclare package Medium = 
          Medium, m_flow(min=if allowFlowReversal and not static then -inf else 0)) 
        "Fluid inlet port" 
                         annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.Ports.FluidPort_b port_b(
                                                 redeclare package Medium = 
          Medium, m_flow(max=if allowFlowReversal and not static then +inf else 0)) 
        "Fluid outlet port" 
                          annotation (extent=[90,-10; 110,10]);
    Medium.BaseProperties[n] medium(
      each preferredMediumStates=if static then false else true,
      p(start=p_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));
      
     annotation (Diagram, Icon(Rectangle(extent=[-100,40; 100,-40], style(
            color=69,
            gradient=2,
            fillColor=69))),
        Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method, such as <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_thermal\">DistributedPipe_thermal</a> and <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_hydraulic\">DistributedPipe_hydraulic</a>. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). Half a momentum balance is used on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
<p><b>Reducing the number of numerical states</b></p>
<p>The default settings of this model result in <tt>n*(2 + nX - 1)</tt> numerical states with <tt>nX</tt> as the number of substances: <tt>n</tt> pressure states, <tt>n</tt> enthalpy or temperature states and <tt>n*nX-1</tt> mass fraction states. Depending on the simulation task the model efficiency may be increased if the number of numerical states is reduced. The following model options exist:</p>
<ul>
<li><tt><b>static</b></tt> - if true, no numerical states are present, static mass and energy balances</li>
<li><tt><b>singleState_hydraulic</b></tt> - if true, only one pressure state is present, just two momentum balances, algebraic constraints for remaining pressures</li>
<li><tt><b>singleState_thermal</b></tt> - if true, only one enthalpy or temperature state is present, just one dynamic energy balance, algebraic constraints for remaining properties</li>
</ul>
Selecting a fixed medium composition along the entire flow path is currently not possible.
<p>
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
    protected 
    final parameter SI.Pressure[np] dp0={if singleState_hydraulic then dp_start else dp_start/(if i
         > 1 and i < np then n else 2*n) for i in 1:np};
    SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
    /*SI.Density d_a=if use_d_nominal then d_nominal else Medium.density_ph(port_a.p, port_a.h);
  SI.Density d_b=if use_d_nominal then d_nominal else Medium.density_ph(port_b.p, port_b.h);
  Above equations currently produce nonlinear systems and numerical Jacobians*/
    SI.Density d_a=if use_d_nominal then d_nominal else d[1];//approximation
    SI.Density d_b=if use_d_nominal then d_nominal else d[n];//approximation
      
  equation 
    // Boundary conditions
    port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium[1].h);
    port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium[n].h);
    port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium[1].Xi);
    port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium[n].Xi);
    port_a.m_flow = m_flow[1];
    port_b.m_flow = -m_flow[n + 1];
      
    // Distributed flow quantities
    for i in 2:n loop
      H_flow[i] = semiLinear(m_flow[i], medium[i - 1].h, medium[i].h);
      mXi_flow[i, :] = semiLinear(m_flow[i], medium[i - 1].Xi, medium[i].Xi);
      v[i] = m_flow[i]/(medium[i - 1].d + medium[i].d)*2/area;
    end for;
    H_flow[1] = port_a.H_flow;
    H_flow[n + 1] = -port_b.H_flow;
    mXi_flow[1, :] = port_a.mXi_flow;
    mXi_flow[n + 1, :] = -port_b.mXi_flow;
    v[1] = m_flow[1]/d_a/area;
    v[n + 1] = m_flow[n + 1]/d_b/area;
      
    // Total quantities
    for i in 1:n loop
      m[i] =Vi[i]*medium[i].d;
      mXi[i, :] = m[i]*medium[i].Xi;
      U[i] = m[i]*medium[i].u;
    end for;
      
    //Mass and energy balances
    if static then
    //steady state mass and energy balances, no numerical states, no flow reversal possible
      for i in 1:n loop
        0 = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        zeros(Medium.nXi) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        0 = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
    elseif singleState_thermal then
    //dynamic mass balances, one dynamic energy balance, n pressure states (if not singleState_hydraulic), 1 "thermal" (h or T) state
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
      end for;
      der(U[nl]) = H_flow[nl] - H_flow[nl + 1] + Qs_flow[nl];
      for i in 1:nl - 1 loop
        medium[i].h = medium[nl].h;
      end for;
      for i in nl + 1:n loop
        medium[i].h = medium[nl].h;
      end for;
    else
    //dynamic mass and energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
      for i in 1:n loop
        der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
        der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
        der(U[i]) = H_flow[i] - H_flow[i + 1] + Qs_flow[i];
      end for;
    end if;
    for i in 1:n loop
      assert((allowFlowReversal and not static) or (m_flow[i] >= 0), "Flow reversal not allowed in distributed pipe");
    end for;
      
  //Momentum Balance, dp contains contributions from acceleration, gravitational and friction effects
  if singleState_hydraulic then //two momentum balances, one on each side of pressure state
      dp[1] = port_a.p - medium[nl].p;
      dp[2] = medium[nl].p - port_b.p;
      if n == 2 then
        medium[2].p = medium[1].p;
      elseif n > 2 then
        medium[1:nl - 1].p = ones(nl - 1)*medium[nl].p;
        medium[nl + 1:n].p = ones(n - nl)*medium[nl].p;
      end if;
    else
      dp[1]=port_a.p-medium[1].p;
      for i in 2:n loop
        dp[i]=medium[i-1].p-medium[i].p;
      end for;
      dp[np]=medium[n].p-port_b.p;
    end if;
      
  initial equation 
    // Initial conditions
    if not static then
      if initType == Types.Init.NoInit then
      // no initial equations
      elseif initType == Types.Init.SteadyState then
      //steady state initialization
        if singleState_thermal then
          if use_T_start then
          der(medium[nl].T) = 0;
        else
          der(medium[nl].h) = 0;
        end if;
        else
        if use_T_start then
          der(medium.T) = zeros(n);
        else
          der(medium.h) = zeros(n);
        end if;
        end if;
        if not (singleState_hydraulic or Medium.singleState) then
          der(medium.p) = zeros(n);
        elseif singleState_hydraulic then
         der(medium[nl].p) = 0;
        end if;
        for i in 1:n loop
          der(medium[i].Xi) = zeros(Medium.nXi);
        end for;
      elseif initType == Types.Init.InitialValues then
      //Initialization with initial values
        if singleState_thermal then
          if use_T_start then
          medium[nl].T = T_start;
        else
          medium[nl].h = h_start;
        end if;
        else
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        end if;
        if not (singleState_hydraulic or Medium.singleState) then
           medium.p=p_start;
        elseif singleState_hydraulic then
         medium[nl].p=p_start[nl];
        end if;
      elseif initType == Types.Init.SteadyStateHydraulic then
      //Steady state initialization for hydraulic states (p)
        if use_T_start then
          medium.T = ones(n)*T_start;
        else
          medium.h = ones(n)*h_start;
        end if;
        if not (singleState_hydraulic or Medium.singleState) then
          der(medium.p) = zeros(n);
        elseif singleState_hydraulic then
          der(medium[nl].p) = 0;
        end if;
      else
        assert(false, "Unsupported initialization option");
      end if;
    end if;
  end PartialDistributedFlow;
    
      partial model PartialLumpedVolume 
      "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
        import Modelica.Constants;
      extends Modelica_Fluid.Interfaces.Records.PartialInitializationParameters;
      replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
      
      parameter Types.FlowDirection.Temp flowDirection=
                Types.FlowDirection.Unidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Dialog(tab="Advanced"));
      Modelica_Fluid.Interfaces.Ports.FluidPort_a port_a(
                                    redeclare package Medium = Medium,
                         m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
        "Fluid inlet port" annotation (extent=[-112,-10; -92,10]);
      Modelica_Fluid.Interfaces.Ports.FluidPort_b port_b(
                                    redeclare package Medium = Medium,
                         m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
        "Fluid outlet port" annotation (extent=[90,-10; 110,10]);
      Medium.BaseProperties medium(preferredMediumStates=true,
                   p(start=p_start), h(start=h_start),
                   T(start=T_start), Xi(start=X_start[1:Medium.nXi]));
      SI.Energy U "Internal energy of fluid";
      SI.Mass m "Mass of fluid";
      SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
      SI.Volume V_lumped "Volume";
      SI.HeatFlowRate Qs_flow "Heat flow across volume boundaries";
      SI.Power Ws_flow "Work flow across boundaries";
      
      annotation (
       Icon(
          Text(extent=[-144, 178; 146, 116], string="%name"),
          Text(
            extent=[-130, -108; 144, -150],
            style(color=0),
            string="V=%V")), Documentation(info="<html>
</html>"),
        Diagram);
      parameter Boolean allowFlowReversal=
         flowDirection == Types.FlowDirection.Bidirectional 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
         annotation(Evaluate=true, Hide=true);
      
      equation 
      // boundary conditions
      port_a.p = medium.p;
      port_b.p = medium.p;
      port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium.h);
      port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium.h);
      port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium.Xi);
      port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium.Xi);
      
      // Total quantities
      m    = V_lumped*medium.d;
      mXi = m*medium.Xi;
      U    = m*medium.u;
      
      // Mass and energy balance
      der(m)   = port_a.m_flow + port_b.m_flow;
      der(mXi) = port_a.mXi_flow + port_b.mXi_flow;
      der(U)= port_a.H_flow + port_b.H_flow + Qs_flow +Ws_flow;
      
      initial equation 
      // Initial conditions
      if initType == Types.Init.NoInit then
        // no initial equations
      elseif initType == Types.Init.InitialValues then
        if not Medium.singleState then
          medium.p = p_start;
        end if;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      elseif initType == Types.Init.SteadyState then
        if not Medium.singleState then
           der(medium.p) = 0;
        end if;
        der(medium.h) = 0;
        der(medium.Xi) = zeros(Medium.nXi);
      elseif initType == Types.Init.SteadyStateHydraulic then
        if not Medium.singleState then
           der(medium.p) = 0;
        end if;
        if use_T_start then
          medium.T = T_start;
        else
          medium.h = h_start;
        end if;
        medium.Xi = X_start[1:Medium.nXi];
      else
        assert(false, "Unsupported initialization option");
      end if;
      end PartialLumpedVolume;
  end ControlVolumes;
  
  package Records 
    partial model PartialInitializationParameters 
      "Define parameter menu to initialize medium in component that has one medium model" 
      import Modelica_Fluid.Types;
      
      replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
      parameter Types.Init.Temp initType=
                Types.Init.NoInit "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_start = Medium.p_default 
        "Start value of pressure" 
        annotation(Dialog(tab = "Initialization"));
      parameter Boolean use_T_start = true 
        "= true, use T_start, otherwise h_start" 
        annotation(Dialog(tab = "Initialization"), Evaluate=true);
      parameter Medium.Temperature T_start=
        if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
        "Start value of temperature" 
        annotation(Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=
        if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
        "Start value of specific enthalpy" 
        annotation(Dialog(tab = "Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Start value of mass fractions m_i/m" 
        annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    end PartialInitializationParameters;
    
    partial model PartialGuessValueParameters 
      "Define parameter menu to initialize guess values of medium in component that has one medium model" 
      import Modelica_Fluid.Types;
      
      replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
      parameter Medium.AbsolutePressure p_start = Medium.p_default 
        "Guess value of pressure" 
        annotation(Dialog(tab = "Guess Value Initialization"));
      parameter Boolean use_T_start = true 
        "= true, use T_start, otherwise h_start" 
        annotation(Dialog(tab = "Guess Value Initialization"), Evaluate=true);
      parameter Medium.Temperature T_start=
        if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
        "Guess value of temperature" 
        annotation(Dialog(tab = "Guess Value Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=
        if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
        "Guess value of specific enthalpy" 
        annotation(Dialog(tab = "Guess Value Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Guess value of mass fractions m_i/m" 
        annotation (Dialog(tab="Guess Value Initialization", enable=Medium.nXi > 0));
    end PartialGuessValueParameters;
  end Records;

  package Valves 
    partial model PartialValve "Base model for valves" 
      import Modelica_Fluid.Types.CvTypes;
      
    parameter Medium.AbsolutePressure pin_start = p_nom 
        "Start value of inlet pressure" 
      annotation(Dialog(tab = "Initialization"));
      
    extends Modelica_Fluid.Interfaces.ControlVolumes.PartialTwoPortTransport(
      medium_a(p(start=pin_start), T(start=T_start),
               h(start=h_start),   Xi(start=X_start[1:Medium.nXi])),
      medium_b(p(start=pout_start), T(start=T_start),
               h(start=h_start),   Xi(start=X_start[1:Medium.nXi])));
      
    parameter CvTypes.Temp CvData = CvTypes.Av "Selection of flow coefficient" 
       annotation(Dialog(group = "Flow Coefficient"));
    parameter SI.Area Av(fixed = if CvData==CvTypes.Av then true else false,
                         start = m_flow_nom/(sqrt(d_nom*dp_nom))*
                                             flowCharacteristic(stemPosition_nom)) = 0 
        "Av (metric) flow coefficient" 
       annotation(Dialog(group = "Flow Coefficient",
                         enable = (CvData==CvTypes.Av)));
    parameter Real Kv(unit="m3/h")=0 "Kv (metric) flow coefficient" 
      annotation(Dialog(group = "Flow Coefficient",
                        enable = (CvData==CvTypes.Kv)));
    parameter Real Cv(unit="USG/min")=0 "Cv (US) flow coefficient" 
      annotation(Dialog(group = "Flow Coefficient",
                        enable = (CvData==CvTypes.Cv)));
    parameter Medium.AbsolutePressure p_nom "Nominal inlet pressure" 
      annotation(Dialog(group="Nominal operating point"));
    parameter SI.Pressure dp_nom "Nominal pressure drop" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Medium.MassFlowRate m_flow_nom "Nominal mass flowrate" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Medium.Density d_nom = 1000 "Nominal inlet density" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Real stemPosition_nom = 1 "Nominal stem position" 
      annotation(Dialog(group="Nominal operating point"));
    parameter Boolean CheckValve=false "Reverse flow stopped";
      
    replaceable function flowCharacteristic = 
        Modelica_Fluid.SubClasses.Valves.ValveCharacteristics.linear 
      extends Modelica_Fluid.Interfaces.Valves.ValveCharacteristics.baseFun 
        "Inherent flow characteristic" 
      annotation(choicesAllMatching=true);
    parameter Medium.AbsolutePressure pout_start = p_nom-dp_nom 
        "Start value of outlet pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
        "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(pin_start,h_start,X_start) 
        "Start value of inlet temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, X_start[1:Medium.nXi]) else Medium.h_default 
        "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
      
    parameter Real delta=0.01 "Regularisation factor" annotation(Dialog(tab="Advanced"));
      
    Modelica.Blocks.Interfaces.RealInput stemPosition 
        "Stem position in the range 0-1" 
                                       annotation (extent=[-20,70; 20,110],   rotation=-90);
      
    Medium.Density d "Density at port a";
    Medium.Temperature T "Temperature at port a";
    protected 
    function sqrtR = Utilities.regRoot(delta = delta*dp_nom);
    annotation (
      Icon(Text(extent=[-143,-66; 148,-106],  string="%name"),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1))),
      Diagram,
      Documentation(info="<HTML>
<p>This is the base model for the <tt>ValveIncompressible</tt>, <tt>ValveVaporizing</tt>, and <tt>ValveCompressible</tt> valve models. The model is based on the IEC 534 / ISA S.75 standards for valve sizing.
<p>The model optionally supports reverse flow conditions (assuming symmetrical behaviour) or check valve operation, and has been suitably modified to avoid numerical singularities at zero pressure drop. 
<p><b>Modelling options</b></p>
<p>The following options are available to specify the valve flow coefficient in fully open conditions:
<ul><li><tt>CvData = Modelica_Fluid.Types.CvTypes.Av</tt>: the flow coefficient is given by the metric <tt>Av</tt> coefficient (m^2).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Kv</tt>: the flow coefficient is given by the metric <tt>Kv</tt> coefficient (m^3/h).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Cv</tt>: the flow coefficient is given by the US <tt>Cv</tt> coefficient (USG/min).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.OpPoint</tt>: the flow is computed from the nominal operating point specified by <tt>p_nom</tt>, <tt>dp_nom</tt>, <tt>m_flow_nom</tt>, <tt>d_nom</tt>, <tt>stemPosition_nom</tt>.
</ul>
<p>The nominal pressure drop <tt>dp_nom</tt> must always be specified; to avoid numerical singularities, the flow characteristic is modified for pressure drops less than <tt>b*dp_nom</tt> (the default value is 1% of the nominal pressure drop). Increase this parameter if numerical problems occur in valves with very low pressure drops.
<p>If <tt>CheckValve</tt> is true, then the flow is stopped when the outlet pressure is higher than the inlet pressure; otherwise, reverse flow takes place.
<p>The inherent flow characteristic <tt>flowCharacteristic</tt>, linear by default, can be replaced by any user-defined function (e.g. equal percentage, quick opening, etc.).
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
      Coordsys(grid=[2,2], scale=0));
    initial equation 
    if CvData == CvTypes.Kv then
      Av = 2.7778e-5*Kv "Unit conversion";
    elseif CvData == CvTypes.Cv then
      Av = 2.4027e-5*Cv "Unit conversion";
    end if;
    assert(CvData>=0 and CvData<=3, "Invalid CvData");
    equation 
    T = medium_a.T;
    d = medium_a.d;
    end PartialValve;

    package ValveCharacteristics 
      partial function baseFun "Base class for valve characteristics" 
        extends Modelica.Icons.Function;
        input Real pos "Stem position (per unit)";
        output Real rc "Relative flow coefficient (per unit)";
      end baseFun;
    end ValveCharacteristics;
  end Valves;

  package FlowMachines 
  partial model PartialPump "Base model for centrifugal pumps" 
      import Modelica.SIunits.Conversions.NonSIunits.*;
      import Modelica.Constants;
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model" 
                     annotation(choicesAllMatching=true);
    Medium.BaseProperties fluid(p(start=pin_start),h(start=h_start)) 
        "Fluid properties at the inlet";
    replaceable package SatMedium = 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium 
        "Saturated medium model (required only for NPSH computation)";
    replaceable function flowCharacteristic = 
        Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.baseFlow 
        "Head vs. q_flow characteristic at nominal speed and density" 
      annotation(Dialog(group="Characteristics"), choicesAllMatching=true);
    parameter Boolean usePowerCharacteristic = false 
        "Use powerCharacteristic (vs. efficiencyCharacteristic)" 
       annotation(Dialog(group="Characteristics"));
    replaceable function powerCharacteristic = 
      Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.basePower 
        "Power consumption vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics", enable = usePowerCharacteristic),
                 choicesAllMatching=true);
    replaceable function efficiencyCharacteristic = 
      Modelica_Fluid.SubClasses.Flowmachines.PumpCharacteristics.constantEfficiency
          (eta_nom=0.8) extends 
        Modelica_Fluid.Interfaces.FlowMachines.PumpCharacteristics.baseEfficiency
        "Efficiency vs. q_flow at nominal speed and density" 
      annotation(Dialog(group="Characteristics",enable = not usePowerCharacteristic),
                 choicesAllMatching=true);
    parameter AngularVelocity_rpm N_nom = 1500 "Nominal rotational speed" 
      annotation(Dialog(group="Characteristics"));
    parameter Medium.Density d_nom = 1000 "Nominal fluid density" 
      annotation(Dialog(group="Characteristics"));
    parameter Integer Np_nom(min=1) = 1 "Nominal number of pumps in parallel";
    parameter SI.Mass M = 0 "Fluid mass inside the pump";
    parameter Boolean checkValve=true "Reverse flow stopped";
    parameter Types.FlowDirection.Temp flowDirection=
                     Types.FlowDirection.Unidirectional 
        "Unidirectional (inlet -> outlet) or bidirectional flow component" 
       annotation(Dialog(tab="Advanced"));
    parameter Boolean computeNPSHa=false "Compute NPSH Available at the inlet";
    parameter Medium.AbsolutePressure pin_start 
        "Guess value for inlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure pout_start 
        "Guess value for outlet pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start = true 
        "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(pin_start,h_start,X_start) 
        "Guess value for temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, X_start) else Medium.h_default 
        "Guess value for specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Guess value for mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter SI.MassFlowRate m_flow_start = 0 
        "Guess value for mass flow rate (total)" 
      annotation(Dialog(tab="Initialization"));
    constant SI.Acceleration g=Modelica.Constants.g_n;
  //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
  //    "Initialisation option";
    Modelica_Fluid.Interfaces.Ports.FluidPort_a inlet(
                                 redeclare package Medium = Medium,
        p(start=pin_start),
        m_flow(start = m_flow_start,
               min = if allowFlowReversal and not checkValve then -Constants.inf else 0)) 
    annotation (extent=[-100,-40; -60,0]);
    Modelica_Fluid.Interfaces.Ports.FluidPort_b outlet(
                                  redeclare package Medium = Medium,
        p(start=pout_start),
        m_flow(start = -m_flow_start,
               max = if allowFlowReversal and not checkValve then +Constants.inf else 0)) 
    annotation (extent=[40,12; 80,52]);
    SI.Pressure dp = outlet.p - inlet.p "Pressure increase";
    SI.Height head = dp/(d*g) "Pump head";
    Medium.Density d "Liquid density at the inlet";
    Medium.SpecificEnthalpy h_out(start=h_start) 
        "Enthalpy of the liquid flowing out of the pump";
    Medium.Temperature Tin "Liquid inlet temperature";
    SI.MassFlowRate m_flow = inlet.m_flow "Mass flow rate (total)";
    SI.MassFlowRate m_flow_single = m_flow/Np "Mass flow rate (single pump)";
    SI.VolumeFlowRate q_flow = m_flow/d "Volume flow rate (total)";
    SI.VolumeFlowRate q_flow_single = q_flow/Np 
        "Volume flow rate (single pump)";
    AngularVelocity_rpm N "Shaft rotational speed";
    Integer Np(min=1) "Number of pumps in parallel";
    SI.Power W_single "Power Consumption (single pump)";
    SI.Power W_tot = W_single*Np "Power Consumption (total)";
    constant SI.Power W_eps=1e-8 
        "Small coefficient to avoid numerical singularities in efficiency computations";
    Real eta "Global Efficiency";
    SI.Length NPSHa "Net Positive Suction Head available";
    Medium.AbsolutePressure pv "Saturation pressure of inlet liquid";
    Real s(start = m_flow_start) 
        "Curvilinear abscissa for the flow curve in parametric form";
    Modelica.Blocks.Interfaces.IntegerInput in_Np 
      annotation (extent=[16,34; 36,54], rotation=-90);
  //  outer Modelica_Fluid.Components.FluidOptions fluidOptions 
  //    "Global default options";
    protected 
   parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
  equation 
    // Number of pumps in parallel
    Np = in_Np;
    if cardinality(in_Np)==0 then
      in_Np = Np_nom "Number of pumps selected by parameter";
    end if;
      
    // Flow equations
    if noEvent(s > 0 or (not checkValve)) then
      // Flow characteristics when check valve is open
      q_flow_single = s;
      // head = (N/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6)));
      head = noEvent((((if abs(N) > 1e-6 then N else 1e-6))/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/((if abs(N) > 1e-6 then N else 1e-6))));
    else
      // Flow characteristics when check valve is closed
      head = (N/N_nom)^2*flowCharacteristic(0) - s;
      q_flow_single = 0;
    end if;
      
    // Power consumption  
    if usePowerCharacteristic then
      W_single = (N/N_nom)^3*(d/d_nom)*powerCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6))) 
          "Power consumption (single pump)";
      eta = (dp*q_flow_single)/(W_single + W_eps) "Hydraulic efficiency";
    else
      eta = efficiencyCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-10)));
      W_single = dp*q_flow/eta;
    end if;
    // Fluid properties
    fluid.p = inlet.p;
    fluid.h = inlet.h;
    fluid.Xi = inlet.Xi;
    d = fluid.d;
    Tin = fluid.T;
      
    // Mass and energy balances
    inlet.m_flow + outlet.m_flow = 0 "Mass balance";
    inlet.mXi_flow + outlet.mXi_flow = zeros(Medium.nXi) 
        "Substance mass balance";
    inlet.H_flow=semiLinear(inlet.m_flow,inlet.h,h_out) 
        "Enthalpy flow at the inlet";
    outlet.H_flow=semiLinear(outlet.m_flow,outlet.h,h_out) 
        "Enthalpy flow at the outlet";
    if M > 0 then
      M * der(h_out) = m_flow_single*(inlet.h - outlet.h) + W_single 
          "Dynamic energy balance (density variations neglected)";
    else
      inlet.H_flow + outlet.H_flow + W_single*Np = 0 "Static energy balance";
    end if;
      
    // NPSH computations
    if computeNPSHa then
      pv=SatMedium.saturationPressure(fluid.T);
      NPSHa=(inlet.p-pv)/(d*Modelica.Constants.g_n);
    else
      pv=0;
      NPSHa=0;
    end if;
  /*
initial equation 
  if initOpt == Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == Choices.Init.Options.steadyState then
    if ThermalCapacity then
      der(h)=0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;
*/
    annotation (
      Icon(
        Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64],
            style(pattern=0, fillColor=74)),
        Ellipse(extent=[-60,40; 60,-80],   style(gradient=3)),
        Polygon(points=[-30,12; -30,-48; 48,-20; -30,12],   style(
            pattern=0,
            gradient=2,
            fillColor=7)),
        Text(extent=[-100,-110; 100,-136], string="%name"),
        Text(extent=[-10,60; 18,40],  string="Np")),
      Diagram,
      Documentation(info="<HTML>
<p>This is the base model for the <tt>Pump</tt> and <tt>
PumpMech</tt> pump models.
<p>The model describes a centrifugal pump, or a group of <tt>Np</tt> identical pumps in parallel. The pump model is based on the theory of kinematic similarity: the pump characteristics are given for nominal operating conditions (rotational speed and fluid density), and then adapted to actual operating condition, according to the similarity equations. 
<p><b>Modelling options</b></p>
<p> The nominal hydraulic characteristic (head vs. volume flow rate) is given by the the replaceable function <tt>flowCharacteristic</tt>. 
<p> The pump energy balance can be specified in two alternative ways:
<ul>
<li><tt>usePowerCharacteristic = false</tt> (default option): the replaceable function <tt>efficiencyCharacteristic</tt> (efficiency vs. volume flow rate in nominal conditions) is used to determine the efficiency, and then the power consumption. The default is a constant efficiency of 0.8.
<li><tt>usePowerCharacteristic = true</tt>: the replaceable function <tt>powerCharacteristic</tt> (power consumption vs. volume flow rate in nominal conditions) is used to determine the power consumption, and then the efficiency.
</ul>
<p>
Several functions are provided in the package <tt>PumpCharacteristics</tt> to specify the characteristics as a function of some operating points at nominal conditions.
<p>Depending on the value of the <tt>checkValve</tt> parameter, the model either supports reverse flow conditions, or includes a built-in check valve to avoid flow reversal.
<p>If the <tt>in_Np</tt> input connector is wired, it provides the number of pumps in parallel; otherwise,  <tt>Np_n</tt> parallel pumps are assumed.</p>
<p>It is possible to take into account the heat capacity of the fluid inside the pump by specifying its mass <tt>M</tt> at nominal conditions; this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flow rate. If zero flow rate conditions are always avoided, this dynamic effect can be neglected by leaving the default value <tt>M = 0</tt>, thus avoiding a fast state variable in the model.
<p>If <tt>computeNPSHa = true</tt>, the available net positive suction head is also computed; this requires a two-phase medium model to provide the fluid saturation pressure.
</HTML>",
        revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
      
  end PartialPump;

  package PumpCharacteristics "Functions for pump characteristics" 
      import NonSI = Modelica.SIunits.Conversions.NonSIunits;
      
    partial function baseFlow "Base class for pump flow characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output SI.Height head "Pump head";
    end baseFlow;
      
    partial function basePower 
        "Base class for pump power consumption characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output SI.Power consumption "Power consumption";
    end basePower;
      
    partial function baseEfficiency "Base class for efficiency characteristics" 
      extends Modelica.Icons.Function;
      input SI.VolumeFlowRate q_flow "Volumetric flow rate";
      output Real eta "Efficiency";
    end baseEfficiency;
      
      
      
      
      
      
  end PumpCharacteristics;
  end FlowMachines;

  package FluidStorage   
  end FluidStorage;

  package HeatTransfer 
    partial model PartialPipeHeatTransfer 
      replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(Dialog(tab="No input", enable=false));
      parameter Integer n(min=1)=1 "Number of pipe segments" annotation(Dialog(tab="No input", enable=false));
      SI.HeatFlowRate[n] Q_flow "Heat flow rates";
      parameter SI.Area A_h "Total heat transfer area" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Length d_h "Hydraulic diameter" annotation(Dialog(tab="No input", enable=false));
      parameter SI.Area A_cross "Cross flow area" annotation(Dialog(tab="No input", enable=false));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
        "Thermal port" 
        annotation (extent=[-20,60; 20,80]);
      SI.Temperature[n] T;
    equation 
      
      annotation (Icon(Ellipse(extent=[-60,64; 60,-56], style(
              color=42,
              rgbcolor={127,0,0},
              gradient=3,
              fillColor=1,
              rgbfillColor={232,0,0})), Text(
            extent=[-38,26; 40,-14],
            style(
              color=42,
              rgbcolor={127,0,0},
              gradient=3,
              fillColor=1,
              rgbfillColor={232,0,0},
              fillPattern=7),
            string="%name")), Documentation(info="<html>
Base class for heat transfer models that can be used in model <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_thermal\">DistributedPipe_thermal</a>.
</html>"));
    end PartialPipeHeatTransfer;
  end HeatTransfer;

  package PressureLosses 
    partial package PartialWallFriction 
      "Partial wall friction characteristic (base package of all wall friction characteristics)" 
      
      annotation (Documentation(info="<html>
 
</html>"));
      
    // Constants to be set in subpackages
      constant Boolean use_eta = true 
        "= true, if eta_a/eta_b are used in function, otherwise value is not used";
      constant Boolean use_roughness = true 
        "= true, if roughness is used in function, otherwise value is not used";
      constant Boolean use_dp_small = true 
        "= true, if dp_small is used in function, otherwise value is not used";
      constant Boolean use_m_flow_small = true 
        "= true, if m_flow_small is used in function, otherwise value is not used";
      constant Boolean dp_is_zero = false 
        "= true, if no wall friction is present, i.e., dp = 0 (function massFlowRate_dp() cannot be used)";
      
    // pressure loss characteristic functions
      replaceable partial function massFlowRate_dp 
        "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction" 
        import SI = Modelica.SIunits;
        extends Modelica.Icons.Function;
        
        input SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input SI.DynamicViscosity eta_a 
          "Dynamic viscosity at port_a (dummy if use_eta = false)";
        input SI.DynamicViscosity eta_b 
          "Dynamic viscosity at port_b (dummy if use_eta = false)";
        input SI.Length length "Length of pipe";
        input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
        input SI.Length roughness(min=0) = 2.5e-5 
          "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
        input SI.AbsolutePressure dp_small = 1 
          "Turbulent flow if |dp| >= dp_small (dummy if use_dp_small = false)";
        
        output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
      annotation (Documentation(info="<html>
 
</html>"));
      end massFlowRate_dp;
      
      replaceable partial function pressureLoss_m_flow 
        "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction" 
        import SI = Modelica.SIunits;
        extends Modelica.Icons.Function;
        
        input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
        input SI.Density d_a "Density at port_a";
        input SI.Density d_b "Density at port_b";
        input SI.DynamicViscosity eta_a 
          "Dynamic viscosity at port_a (dummy if use_eta = false)";
        input SI.DynamicViscosity eta_b 
          "Dynamic viscosity at port_b (dummy if use_eta = false)";
        input SI.Length length "Length of pipe";
        input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
        input SI.Length roughness(min=0) = 2.5e-5 
          "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
        input SI.MassFlowRate m_flow_small = 0.01 
          "Turbulent flow if |m_flow| >= m_flow_small (dummy if use_m_flow_small = false)";
        output SI.Pressure dp "Pressure drop (dp = port_a.p - port_b.p)";
        
      annotation (Documentation(info="<html>
 
</html>"));
      end pressureLoss_m_flow;
      
    end PartialWallFriction;

    model PartialHydraulicResistance 
      "Generic pressure drop component with constant turbulent loss factor data and without an icon" 
      
      extends Modelica_Fluid.Interfaces.Records.PartialGuessValueParameters;
      extends Modelica_Fluid.Interfaces.ControlVolumes.PartialTwoPortTransport(
                    medium_a(p(start=p_start), h(start=h_start),
                             T(start=T_start), Xi(start=X_start[1:Medium.nXi])),
                    medium_b(p(start=p_start), h(start=h_start),
                             T(start=T_start), Xi(start=X_start[1:Medium.nXi])));
      
      SI.ReynoldsNumber Re = Utilities.ReynoldsNumber_m_flow(
            m_flow, (Medium.dynamicViscosity(medium_a) + Medium.dynamicViscosity(medium_b))/2,
            data.D_Re) if show_Re "Reynolds number at diameter data.D_Re";
      parameter 
        Modelica_Fluid.SubClasses.PressureLosses.QuadraticTurbulent.LossFactorData
        data "Loss factor data";
      parameter Boolean show_Re = false 
        "= true, if Reynolds number is included for plotting" 
         annotation (Evaluate=true, Dialog(tab="Advanced"));
      parameter Boolean from_dp = true 
        "= true, use m_flow = f(dp) else dp = f(m_flow)" 
        annotation (Evaluate=true, Dialog(tab="Advanced"));
      parameter Boolean use_Re = false 
        "= true, if turbulent region is defined by Re, otherwise by dp_small or m_flow_small"
        annotation(Evaluate=true, Dialog(tab="Advanced"));
      parameter SI.AbsolutePressure dp_small = 1 
        "Turbulent flow if |dp| >= dp_small" 
        annotation(Dialog(tab="Advanced", enable=not use_Re and from_dp));
      parameter SI.MassFlowRate m_flow_small = 0.01 
        "Turbulent flow if |m_flow| >= m_flow_small" 
        annotation(Dialog(tab="Advanced", enable=not use_Re and not from_dp));
      
      annotation (
        Diagram,
        Icon,
        Documentation(info="<html>
<p>
This model computes the pressure loss of a pipe
segment (orifice, bending etc.) with a minimum amount of data
provided via parameter <b>data</b>.
If available, data should be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
</p>
<p>
The following equations are used:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
        Re = |v|*D*&rho;/&eta;
</pre>
<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>flow type</b></td>
    <td><b>&zeta;</b> = </td>
    <td><b>flow region</b></td></tr>
<tr><td>turbulent</td>
    <td><b>zeta1</b> = const.</td>
    <td>Re &ge;  Re_turbulent, v &ge; 0</td></tr>
<tr><td></td>
    <td><b>zeta2</b> = const.</td>
    <td>Re &ge; Re_turbulent, v &lt; 0</td></tr>
<tr><td>laminar</td>
    <td><b>c0</b>/Re</td>
    <td>both flow directions, Re small; c0 = const.</td></tr>
</table>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<br>
     When the Reynolds number Re is below \"Re_turbulent\", the
     flow is laminar for small flow velocities. For higher 
     velocities there is a transition region from 
     laminar to turbulent flow. The loss factor for
     laminar flow at small velocities is defined by the often occuring
     approximation c0/Re. If c0 is different for the two
     flow directions, the mean value has to be used 
     (c0 = (c0_ab + c0_ba)/2).<li>
<li> The equation \"&Delta;p = 0.5*&zeta;*&rho;*v*|v|\" is either with
     respect to port_a or to port_b, depending on the definition
     of the particular loss factor &zeta; (in some references loss
     factors are defined with respect to port_a, in other references
     with respect to port_b).</li>
 
<li> Re = |v|*D_Re*&rho;/&eta; = |m_flow|*D_Re/(A_Re*&eta;) 
     is the Reynolds number at the smallest cross
     section area. This is often at port_a or at port_b, but can
     also be between the two ports. In the record, the diameter
     D_Re of this smallest cross section area has to be provided, as
     well, as Re_turbulent, the absolute value of the 
     Reynolds number at which
     the turbulent flow starts. If Re_turbulent is different for
     the two flow directions, use the smaller value as Re_turbulent.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> A is the cross section area with A = &pi;(D/2)^2.
<li> &eta; is the dynamic viscosity.</li>
</ul>
<p>
The laminar and the transition region is usually of
not much technical interest because the operating point is
mostly in the turbulent regime. For simplification and for
numercial reasons, this whole region is described by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. The polynomials start at 
Re = |m_flow|*4/(&pi;*D_Re*&eta;), where D_Re is the
smallest diameter between port_a and port_b.
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"c0/Re\". Note, the pressure drop
equation above in the laminar region is always defined
with respect to the smallest diameter D_Re.
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));
    equation 
      if from_dp then
        m_flow = if use_Re then 
          Modelica_Fluid.SubClasses.PressureLosses.QuadraticTurbulent.massFlowRate_dp_and_Re(
              dp, 
              medium_a.d, 
              medium_b.d, 
              Medium.dynamicViscosity(medium_a.state), 
              Medium.dynamicViscosity(medium_b.state), 
              data) else 
          Modelica_Fluid.SubClasses.PressureLosses.QuadraticTurbulent.massFlowRate_dp(
              dp, 
              medium_a.d, 
              medium_b.d, 
              data, 
              dp_small);
      else
        dp = if use_Re then 
          Modelica_Fluid.SubClasses.PressureLosses.QuadraticTurbulent.pressureLoss_m_flow_and_Re(
              m_flow, 
              medium_a.d, 
              medium_b.d, 
              Medium.dynamicViscosity(medium_a.state), 
              Medium.dynamicViscosity(medium_b.state), 
              data) else 
          Modelica_Fluid.SubClasses.PressureLosses.QuadraticTurbulent.pressureLoss_m_flow(
              m_flow, 
              medium_a.d, 
              medium_b.d, 
              data, 
              m_flow_small);
      end if;
    end PartialHydraulicResistance;
  end PressureLosses;

  package Ports 
    connector FluidPort 
      "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)" 
      
      replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium model" annotation (choicesAllMatching=true);
      
      Medium.AbsolutePressure p "Pressure in the connection point";
      flow Medium.MassFlowRate m_flow 
        "Mass flow rate from the connection point into the component";
      
      Medium.SpecificEnthalpy h 
        "Specific mixture enthalpy in the connection point";
      flow Medium.EnthalpyFlowRate H_flow 
        "Enthalpy flow rate into the component (if m_flow > 0, H_flow = m_flow*h)";
      
      Medium.MassFraction Xi[Medium.nXi] 
        "Independent mixture mass fractions m_i/m in the connection point";
      flow Medium.MassFlowRate mXi_flow[Medium.nXi] 
        "Mass flow rates of the independent substances from the connection point into the component (if m_flow > 0, mXi_flow = m_flow*Xi)";
      
      Medium.ExtraProperty C[Medium.nC] 
        "properties c_i/m in the connection point";
      flow Medium.ExtraPropertyFlowRate mC_flow[Medium.nC] 
        "Flow rates of auxiliary properties from the connection point into the component (if m_flow > 0, mC_flow = m_flow*C)";
      
    end FluidPort;

    connector FluidPort_a "Fluid connector with filled icon" 
      extends Modelica_Fluid.Interfaces.Ports.FluidPort;
      annotation (defaultComponentName="port_a",
                  Diagram(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
                 fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
                 fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")),
           Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
                fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
                fillColor=69))));
    end FluidPort_a;

    connector FluidPort_b "Fluid connector with outlined icon" 
      extends Modelica_Fluid.Interfaces.Ports.FluidPort;
      annotation (defaultComponentName="port_b",
                  Diagram(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
                 fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
                 fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
                 fillColor=7)), Text(extent=[-88, 192; 112, 98], string="%name")),
           Icon(Ellipse(extent=[-100, 100; 100, -100], style(color=69,
                fillColor=69)), Ellipse(extent=[-100, 100; 100, -100], style(color=16,
                fillColor=69)), Ellipse(extent=[-80, 80; 80, -80], style(color=69,
                 fillColor=7))));
    end FluidPort_b;
  end Ports;

  package Sensors 
    partial model PartialAbsoluteSensor 
      "Partial component to model a sensor that measures a potential variable" 
      
      replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium "Medium in the sensor" annotation (
          choicesAllMatching =                                                                        true);
      Modelica_Fluid.Interfaces.Ports.FluidPort_a port(
                                  redeclare package Medium = Medium) 
        annotation (extent=[-10,-110; 10,-90],    rotation=90);
      
      annotation (Documentation(info="<html>
<p>
Partial component to model an <b>absolute sensor</b>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <tt>PartialFlowSensor</tt> instead.
as signal.
</p>
</html>"),
        Diagram,
        Coordsys(grid=[1,1], scale=0));
    equation 
      port.m_flow = 0;
      port.H_flow = 0;
      port.mXi_flow = zeros(Medium.nXi);
    end PartialAbsoluteSensor;

    partial model PartialFlowSensor 
      "Partial component to model sensors that measure flow properties" 
      
      import Modelica.Constants;
      
      replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium "Medium in the sensor"  annotation (
          choicesAllMatching = true);
      Medium.SpecificEnthalpy h "enthalpy in flow";
      Medium.MassFraction[Medium.nXi] Xi "flow composition";
      
      Modelica_Fluid.Interfaces.Ports.FluidPort_a port_a(
                                    redeclare package Medium = Medium,
                         m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
        annotation (extent=[-110,-10; -90,10]);
      Modelica_Fluid.Interfaces.Ports.FluidPort_b port_b(
                                    redeclare package Medium = Medium,
                         m_flow(max=if allowFlowReversal then +Constants.inf else 0)) 
        annotation (extent=[110,-10; 90,10]);
      
      parameter Modelica_Fluid.Types.FlowDirection.Temp flowDirection=
                Modelica_Fluid.Types.FlowDirection.Unidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Dialog(tab="Advanced"));
      
      annotation (Documentation(info="<html>
<p>
Partial component to model a <b>sensor</b> that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.<br>
The model includes zero-volume balance equations. Sensor models inheriting from
this partial class should add a medium instance to calculate the measured property.
</p>
</html>"),
        Diagram,
        Coordsys(grid=[1,1], scale=0));
      
    protected 
      parameter Boolean allowFlowReversal=
         flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "= false, if flow only from port_a to port_b, otherwise reversing flow allowed"
         annotation(Evaluate=true, Hide=true);
    equation 
      port_a.p   = port_b.p;
      // Local *zero-volume* enthalpy and composition
      port_a.H_flow = semiLinear(port_a.m_flow,port_a.h,h);
      port_b.H_flow = semiLinear(port_b.m_flow,port_b.h,h);
      port_a.mXi_flow = semiLinear(port_a.m_flow,port_a.Xi,Xi);
      port_b.mXi_flow = semiLinear(port_b.m_flow,port_b.Xi,Xi);
      // Static balances
      0 = port_a.m_flow + port_b.m_flow;
      0 = port_a.H_flow + port_b.H_flow;
      zeros(Medium.nXi) = port_a.mXi_flow + port_b.mXi_flow;
    end PartialFlowSensor;
  end Sensors;

  package Sources 
  partial model PartialSource 
      "Partial component source with one fluid connector" 
      import Modelica.Constants;
    replaceable package Medium = PackageMedium extends 
        Modelica.Media.Interfaces.PartialMedium 
        "Medium model within the source" 
       annotation (choicesAllMatching=true);
    Modelica_Fluid.Interfaces.Ports.FluidPort_b port(
                                redeclare package Medium = Medium,
                     m_flow(min=if allowFlowReversal then -Constants.inf else 0)) 
      annotation (extent=[90,-10; 110,10],    rotation=0);
    Medium.BaseProperties medium "Medium in the source";
    parameter Types.FlowDirection.Temp flowDirection=
                     Types.FlowDirection.Unidirectional 
        "Unidirectional (out of port_b) or bidirectional flow component" 
                                                                annotation(Dialog(tab="Advanced"));
    protected 
      parameter Boolean allowFlowReversal=
       flowDirection == Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "= false, if flow only out of port_b, otherwise reversing flow allowed"
       annotation(Evaluate=true, Hide=true);
  equation 
    port.p = medium.p;
    port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.Xi);
    annotation (Documentation(info="<html>
<p>
Partial component to model the <b>volume interface</b> of a <b>source</b>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= port.p) is identical to the
     pressure in the volume (= medium.p).</li>
<li> The enthalpy flow rate (= port.H_flow) and the mass flow rates of the
     substances (= port.mX_flow) depend on the direction of the mass flow rate.</li>
</ul>
</html>"),
      Diagram,
      Coordsys(grid=[1,1], scale=0));
  end PartialSource;
  end Sources;
end Interfaces;
