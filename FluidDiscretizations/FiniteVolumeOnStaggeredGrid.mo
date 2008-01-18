package FiniteVolumeOnStaggeredGrid 
  "Discretizations using the finite volume method on a staggered grid" 
  extends Icons.FluidDiscretization;
  
  package BaseClasses "Base classes for pipe models" 
    package FiniteVolumeOnStaggeredGrid 
      "Finite volume method on staggered grid" 
      extends Interfaces.PartialFluidDiscretization;
      redeclare replaceable partial model extends PartialDistributedFlow 
        "Constant cross section finite volume method with momentum balance on a staggered grid" 
        
      // Taken one-to-one from Modelica_Fluid
        
      // Covered herein:
      //   Internal control volumes
      //   Internal flows
      //   Initialization
      //   Geometry
        
        import Modelica_Fluid.Types;
        
        replaceable package Medium = 
            Modelica.Media.Interfaces.PartialMedium "Fluid medium model" 
            annotation (choicesAllMatching=true);
        
      //Discretization
        parameter Integer n(min=1) = 1 "Number of pipe segments";
        
        parameter Boolean use_approxPortProperties=false 
          "=true, port properties for pressure drop correlation are taken from neighboring control volume"
                                                                                                            annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
        
      //Initialization (TODO: Initialization of m_flow for dynamic momentum)
        parameter Types.Init.Temp initType=Types.Init.NoInit 
          "Initialization option" 
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
              p_a_start,
              p_b_start + (p_a_start - p_b_start)/n,
              n) else {(p_a_start + p_b_start)/2} "Start value of pressure";
        parameter Medium.Temperature T_start=if use_T_start then Medium.T_default else 
                  Medium.temperature_phX(
              (p_a_start + p_b_start)/2,
              h_start,
              X_start) "Start value of temperature" 
                annotation(Evaluate=true, Dialog(tab = "Initialization", enable = use_T_start));
        parameter Medium.SpecificEnthalpy h_start=if use_T_start then Medium.specificEnthalpy_pTX(
              (p_a_start + p_b_start)/2,
              T_start,
              X_start) else Medium.h_default "Start value of specific enthalpy"
                annotation(Evaluate=true, Dialog(tab = "Initialization", enable = not use_T_start));
        parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
          "Start value of mass fractions m_i/m" 
              annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
        parameter Medium.MassFlowRate mflow_start 
          "Start value for mass flow rate"                                             annotation(Evaluate=true, Dialog(tab = "Initialization"));
        final parameter SI.Pressure dp_start=p_a_start - p_b_start;
      protected 
        final parameter SI.Pressure[n + 1] dp0=fill(dp_start, n + 1) 
          "pressure difference start values";
        
      //Geometry
      public 
        parameter Boolean isCircular=true 
          "= true if cross sectional area is circular" 
                annotation (Evaluate, Dialog(tab="General", group="Geometry"));
        parameter SI.Length length "Length"       annotation(Dialog(tab="General", group="Geometry"));
        parameter SI.Diameter diameter "Diameter of circular pipe"          annotation(Dialog(group="Geometry", enable=isCircular));
        parameter SI.Length perimeter=if isCircular then Modelica.Constants.pi*diameter else 0 
          "Inner perimeter"                                                                             annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
        parameter SI.Area area=if isCircular then Modelica.Constants.pi*diameter*diameter/4 else 0 
          "Inner cross section area"                annotation(Dialog(tab="General", group="Geometry", enable=not isCircular));
        SI.Volume[n] Vi "Discretized volume, determine in inheriting class ";
        
      // Density approximation
        parameter Boolean use_d_nominal=false 
          "= true, if d_nominal is used, otherwise computed from medium"                                    annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
        parameter SI.Density d_nominal 
          "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)"    annotation(Dialog(tab="Advanced", group="Momentum balance",enable=use_d_nominal));
        SI.Pressure[n + 1] dp(start=dp0) 
          "pressure difference across staggered grid (1 and n+1 between port pressures and medium pressures, others between CVs)";
        
      //Total quantities
        SI.Energy[n] U "Internal energy of fluid";
        SI.Mass[n] m "Fluid mass";
        SI.Mass[n,Medium.nXi] mXi "Substance mass";
        
      //Flow quantities
        inner Medium.MassFlowRate[n + 1] m_flow(each start=mflow_start, each fixed
            = false) "Mass flow rates of fluid across segment boundaries";
        SI.Velocity[n] v 
          "Velocity at CV centroids (used for the stagnation properties)";
      protected 
        SI.Velocity[n+1] v_staggered 
          "Velocity at segment boundaries (used for the stagnation properties)";
      public 
        Medium.MassFlowRate[n + 1,Medium.nXi] mXi_flow 
          "Independent substance mass flow rates across segment boundaries";
        Medium.EnthalpyFlowRate[n + 1] H_flow 
          "Enthalpy flow rates of fluid across segment boundaries";
        
        Medium.BaseProperties[n] medium(
          each preferredMediumStates=true,
          p(start=p_start),
          each h(start=h_start),
          each T(start=T_start),
          each Xi(start=X_start[1:Medium.nXi]));
        
      //Source terms, have to be set in inheriting class (to zero if not used)
      protected 
        Medium.MassFlowRate[n] ms_flow "Mass flow rate, source or sink";
        Medium.MassFlowRate[n,Medium.nXi] msXi_flow 
          "Independent mass flow rates, source or sink";
        SI.HeatFlowRate[n] Qs_flow "Heat flow rate, source or sink";
        SI.Power[n] Ws_flow "Mechanical power, p*der(V) etc.";
        
        SI.Density[n] d=if use_d_nominal then ones(n)*d_nominal else medium.d;
        
        Medium.AbsolutePressure port_a_p "Interface to port pressure";
        Medium.AbsolutePressure port_b_p "Interface to port pressure";
        
      public 
        Medium.DynamicViscosity eta_a "Interface to port dynamic viscosity";
        Medium.DynamicViscosity eta_b "Interface to port dynamic viscosity";
        Medium.Density d_a "Interface to port density";
        Medium.Density d_b "Interface to port density";
        
        parameter Boolean useStagantionProperties = false 
          "Use of stagnation properties in the energy balance 1?"                                                 annotation(Evaluate=true);
        
      equation 
      // Distributed flow quantities, upwind discretization
        for i in 2:n loop
          H_flow[i] = noEvent(if m_flow[i] > 0 then m_flow[i]*medium[i - 1].h else 
            m_flow[i]*medium[i].h);
          mXi_flow[i, :] = noEvent(if m_flow[i] > 0 then m_flow[i]*medium[i - 1].Xi else 
                  m_flow[i]*medium[i].Xi);
        end for;
        for i in 1:n loop
          v[i] = (m_flow[i]+m_flow[i+1])/2/(medium[i].d*area);
        end for;
        
      // Total quantities
        for i in 1:n loop
          m[i] = Vi[i]*medium[i].d;
          mXi[i, :] = m[i]*medium[i].Xi;
          U[i] = m[i]*(medium[i].u + (if useStagantionProperties then v[i]^2/2 else 0));
        end for;
        
      //Mass and energy balances
      //dynamic mass and energy balances, n "thermal" states, n pressure states (if not singleState_hydraulic)
        for i in 1:n loop
          der(m[i]) = m_flow[i] - m_flow[i + 1] + ms_flow[i];
          der(mXi[i, :]) = mXi_flow[i, :] - mXi_flow[i + 1, :] + msXi_flow[i, :];
          der(U[i]) = (H_flow[i] + (if useStagantionProperties then v_staggered[i]^2/2 else 0)) - (H_flow[i + 1] + (if useStagantionProperties then v_staggered[i+1]^2/2 else 0)) + Qs_flow[i];
        end for;
        
      //Pressure difference between control volumes
        dp[1] = port_a_p - medium[1].p;
        for i in 2:n loop
          dp[i] = medium[i - 1].p - medium[i].p "Between control volumes";
        end for;
        dp[n + 1] = medium[n].p - port_b_p;
        
        // flow velocities at boundaries only for stagnation properties
        v_staggered[1] = if useStagantionProperties then m_flow[1]/(d_a*area) else 0;
        for i in 2:n loop
          v_staggered[i] = if useStagantionProperties then m_flow[i]/((medium[i - 1].d + medium[i].d)/2*area) else 0;
        end for;
        v_staggered[n + 1] = if useStagantionProperties then m_flow[n + 1]/(d_b * area) else 0;
        
      initial equation 
      // Initial conditions
        if initType == Types.Init.NoInit then
      // no initial equations
        elseif initType == Types.Init.SteadyState then
        //steady state initialization
          if use_T_start then
            der(medium.T) = zeros(n);
          else
            der(medium.h) = zeros(n);
          end if;
          if not (Medium.singleState) then
            der(medium.p) = zeros(n);
          end if;
          for i in 1:n loop
            der(medium[i].Xi) = zeros(Medium.nXi);
          end for;
        elseif initType == Types.Init.InitialValues then
        //Initialization with initial values
          if use_T_start then
            medium.T = ones(n)*T_start;
          else
            medium.h = ones(n)*h_start;
          end if;
          if not Medium.singleState then
            medium.p = p_start;
          end if;
        elseif initType == Types.Init.SteadyStateHydraulic then
        //Steady state initialization for hydraulic states (p)
          if use_T_start then
            medium.T = ones(n)*T_start;
          else
            medium.h = ones(n)*h_start;
          end if;
          if not Medium.singleState then
            der(medium.p) = zeros(n);
          end if;
        else
          assert(false, "Unsupported initialization option");
        end if;
        
        annotation (
          Diagram,
          Icon(Rectangle(extent=[-100,40; 100,-40], style(
                color=69,
                gradient=2,
                fillColor=69))),
          Documentation(info="<html>
<p>The model <b>PartialDistributedFlow</b> is used as a base class for pipe flows with one-dimensional spatial discretization according to the finite volume method. The flow path is divided into <tt><b>n</b></tt> segments.</p>
 
<p><b>Mass and energy balances</b></p>
<p>One total mass and one energy balance is formed across each segment. If the medium contains more than one component, substance mass balances are added. Changes in potential and kinetic energy are neglected in the energy balance. The following source (or sink) terms are used in the balances and must be specified in extending models to complete this partial class:</p>
<ul>
<li>Energy balance: <tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and <tt><b>Ws_flow</b></tt>, e.g. mechanical power</li>
<li>Total mass balance: <tt><b>ms_flow</b></tt>, e.g. condensing mass flow of negligible volume such as water in moist air</li>
<li>Substance mass balance: <tt><b>msXi_flow</b></tt>, as above</li>
</ul>
If the flag <tt>static</tt> is <b>true</b> then no mass or energy is stored in the component and the mass and energy balances are reduced to a quasi steady-state formulation. It should be noted that dynamic balances are required if flow reversal should be allowed.
<p>In addition the volume vector <tt><b>Vi</b></tt>, which specifies the volume of each segment and the pressure drop (or rise) in each segment <tt><b>dp</b></tt> must be provided in the extending class.
 
<p><b>Momentum balance</b></p>
<p>The momentum balance is always static, i.e. no dynamic momentum term is used. The momentum balances are formed across the segment boundaries (staggered grid). The default symmetric model is characterized by half a momentum balance on each end of the flow model resulting in a total of n-1 full and 2 half momentum balances. Connecting two pipes therefore results in an algebraic pressure at the ports. Specifying a good start value for the port pressure is essential in order to solve large systems. Non-symmetric variations are obtained by chosing a different value for the parameter <tt><b>modelStructure</b></tt>. Options include:
<ul>
<li><tt>a_v_b</tt>: default setting with two half momentum balances</li>
<li><tt>av_b</tt>: full momentum balance between nth volume and <tt>port_b</tt>, potential pressure state at <tt>port_a</tt></li>
<li><tt>a_vb</tt>: full momentum balance between first volume and <tt>port_a</tt>, potential pressure state at <tt>port_b</tt></li>
<li><tt>avb</tt>: n-1 momentum balances between first and nth volume, potential pressure states at both ports. It's use should be avoided, since not the entire pipe length is taken into account.
</ul></p>
 
<p>The term <tt>dp</tt> is unspecified in this partial class. When extending from this model it may contain
<ul>
<li>pressure drop due to friction and other dissipative losses</li>
<li>changes in pressure resulting from significant variation of flow velocity along the flow path (with the assumption of a constant cross sectional area it must result from fluid density changes, such as in two-phase flow)</li>
<li>gravity effects for non-horizontal pipes</li>
</ul>
At least one relationship between pressure difference and massflow rate (dp=dp(m_flow)) is required in the extending class for this term to receive a fully determined model.
 
When connecting two components, e.g. two pipes, the momentum balance across the connection point reduces to</p> 
<pre>pipe1.port_b.p = pipe2.port_a.p</pre>
<p>This is only true if the flow velocity remains the same on each side of the connection. For any significant change in diameter (and if the resulting effects, such as change in kinetic energy, cannot be neglected) an adapter component should be used. This also allows for taking into account friction losses with respect to the actual geometry of the connection point.</p>
 
 
 
 
 
</html>",     revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
        
      end PartialDistributedFlow;
    end FiniteVolumeOnStaggeredGrid;
    extends Icons.BaseClassLibrary;
    
    package FiniteVolumeOnStaggeredGridWithCorrelations 
      "Adds pressure drop and heat transfer correlations" 
      extends FiniteVolumeOnStaggeredGrid;
      redeclare replaceable partial model extends PartialDistributedFlow 
      // Taken one-to-one from Modelica_Fluid
        
      // Covered herein:
      //   Pressure loss correlation
      //   Heat transfer correlation
      //   Internal friction pressure losses
      //   Heat transfer
        
        replaceable package Medium = 
            Modelica.Media.Interfaces.PartialMedium "Fluid medium model" 
            annotation (choicesAllMatching=true);
        
        parameter SI.Area area_h=perimeter*length "Heat transfer area" 
                                                           annotation(Dialog(tab="General", group="Heat transfer"));
        final parameter SI.Volume V=area*length;
        parameter SI.Length height_ab=0.0 "Height(port_b) - Height(port_a)" 
                                                                 annotation(Dialog(tab="General", group="Geometry"),Evaluate=true);
        
      //Pressure Drop
        replaceable package WallFriction = 
            PressureLosses.WallFrictionCorrelations.QuadraticTurbulent extends 
          PressureLosses.WallFrictionCorrelations.PartialWallFriction 
          "Characteristic of wall friction"  annotation(Dialog(tab="General", group="Pressure loss"),choicesAllMatching=true);
        parameter Boolean use_d_nominal=false 
          "= true, if d_nominal is used, otherwise computed from medium"                        annotation(Dialog(tab="Advanced", group="Momentum balance"),Evaluate=true);
        parameter SI.Diameter d_h=4*area/perimeter "Hydraulic diameter" 
                                                               annotation(Dialog(tab="General", group="Pressure loss"));
        parameter SI.Length roughness(min=0) = 2.5e-5 
          "Absolute roughness of pipe (default = smooth steel pipe)" 
                                                                   annotation (
          Dialog(
          tab="General",
          group="Pressure loss",
          enable=WallFriction.use_roughness));
        parameter Boolean use_eta_nominal=false 
          "= true, if eta_nominal is used, otherwise computed from medium"                    annotation(Dialog(tab="General", group="Pressure loss"),Evaluate=true);
        parameter SI.DynamicViscosity eta_nominal 
          "Nominal dynamic viscosity (e.g. eta_liquidWater = 1e-3, eta_air = 1.8e-5)"
                                                                        annotation(Dialog(tab="General", group="Pressure loss",enable=use_nominal));
        parameter Boolean show_Re=false 
          "= true, if Reynolds number is included for plotting" 
      annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
        parameter Boolean from_dp=true 
          " = true, use m_flow = f(dp), otherwise dp = f(m_flow)" 
      annotation (Evaluate=true, Dialog(tab="Advanced", group="Pressure loss"));
        parameter SI.AbsolutePressure dp_small=1 
          "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
                                                                                                          annotation (Dialog(
        tab="Advanced",
        group="Pressure loss",
        enable=from_dp and WallFriction.use_dp_small));
        
        parameter SI.MassFlowRate m_flow_small=0.01 
          "Turbulent flow if |m_flow| >= m_flow_small (only used if WallFriction=QuadraticTurbulent)"
                                                                                                          annotation (
        Dialog(
        tab="Advanced",
        group="Pressure loss",
        enable=not from_dp and WallFriction.use_m_flow_small));
        
        SI.ReynoldsNumber[n] Re=Modelica_Fluid.Utilities.ReynoldsNumber_m_flow(
            m_flow,
            (eta_a + eta_b)/2,
            diameter) if                                                                                     show_Re 
          "Reynolds number of pipe flow";
        inner Medium.ThermodynamicState[n] state;
        
      // Heat transfer correlation
        replaceable 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PipeHT_constAlpha 
          heatTransfer(
          redeclare final package Medium = Medium,
          final n=n,
          final d_h=d_h,
          final A_h=area_h,
          final A_cross=area,
          final L=length,
          T=medium.T) extends 
          Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer(
          redeclare final package Medium = Medium,
          final n=n,
          final d_h=d_h,
          final A_h=area_h,
          final A_cross=area,
          final L=length,
          T=medium.T) "Convective heat transfer" 
          annotation (Dialog(tab="General", group="Heat transfer"),editButton=true,choicesAllMatching, extent=[-20,-20;20,20]);
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[n] thermalPort 
          "Thermal port" 
                   annotation (extent=[-10,44; 10,64]);
        outer Modelica_Fluid.Ambient ambient "Ambient conditions";
        
      protected 
        SI.Pressure[n + 1] dp_friction;
        
        SI.DynamicViscosity[n] eta=if not WallFriction.use_eta then fill(1.e-10, n) else 
                  (if use_eta_nominal then fill(eta_nominal, n) else 
            Medium.dynamicViscosity(medium.state));
        
      equation 
        // for heat transfer
        state = medium.state;
        
        //Pressure drop and gravity
        if from_dp and not WallFriction.dp_is_zero then
          for i in 2:n loop
            m_flow[i] = WallFriction.massFlowRate_dp(
              dp_friction[i] - height_ab/n*ambient.g*(d[i - 1] + d[i])/2,
              d[i - 1],
              d[i],
              eta[i - 1],
              eta[i],
              length/n,
              d_h,
              roughness,
              dp_small);
          end for;
        else
          for i in 2:n loop
            dp_friction[i] = WallFriction.pressureLoss_m_flow(
              m_flow[i],
              d[i - 1],
              d[i],
              eta[i - 1],
              eta[i],
              length/n,
              d_h,
              roughness,
              m_flow_small) + height_ab/n*ambient.g*(d[i - 1] + d[i])/2;
          end for;
        end if;
        connect(thermalPort, heatTransfer.thermalPort) 
        annotation (points=[0,54; 0,14], style(color=42, rgbcolor={191,0,0}));
        annotation (
          Icon(
            Rectangle(extent=[-100,44; 100,40], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=10,
                rgbfillColor={95,95,95})),
            Rectangle(extent=[-100,-40; 100,-44], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=10,
                rgbfillColor={95,95,95})),
            Ellipse(extent=[-82,10; -62,-10], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=0,
                rgbfillColor={0,0,0})),
            Ellipse(extent=[-34,10; -14,-10], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=0,
                rgbfillColor={0,0,0})),
            Ellipse(extent=[14,10; 34,-10], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=0,
                rgbfillColor={0,0,0})),
            Ellipse(extent=[62,10; 82,-10], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=0,
                rgbfillColor={0,0,0})),
            Text(
              extent=[-143,-42; 157,-92],
              string="%name",
              style(gradient=2, fillColor=69))),
          Diagram,
          Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
 
</html>",       revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
        
      end PartialDistributedFlow;
    end FiniteVolumeOnStaggeredGridWithCorrelations;
  end BaseClasses;
  
 package StaticMomentum "Discretizations with static momentum balances" 
  extends Icons.VariantLibrary;
   package Symmetric 
      "Modelica_Fluid inspired discretization using the finite volume method on a staggered grid with static momentum balance (symmetric)" 
         //extends Interfaces.PartialFluidDiscretization(isSymmetric=true);
      
     extends BaseClasses.FiniteVolumeOnStaggeredGridWithCorrelations(isSymmetric=true);
      
     redeclare replaceable partial model extends PartialDistributedFlow 
        
       // Taken one-to-one from Modelica_Fluid
        
       // Covered herein:
       //   TBD
        
     equation 
       //Static momentum balance
       0 = - (medium[1].p - port_a_p) - dp_friction[1];
       for i in 2:n loop
         0 = - (medium[i].p - medium[i - 1].p) - dp_friction[i];
       end for;
       0 = - (port_b_p - medium[n].p) - dp_friction[n + 1];
        
       // Source terms  
       Qs_flow = heatTransfer.Q_flow;
       Ws_flow = zeros(n);
       ms_flow = zeros(n);
       msXi_flow = zeros(n, Medium.nXi);
       Vi = ones(n)*V/n;
       //Pressure drop and gravity
       if from_dp and not WallFriction.dp_is_zero then
         m_flow[1] = WallFriction.massFlowRate_dp(
           dp_friction[1] - height_ab/2/n*ambient.g*d[1],
           d_a,
           d[1],
           eta_a,
           eta[1],
           length/n/2,
           d_h,
           roughness,
           dp_small);
         m_flow[n + 1] = WallFriction.massFlowRate_dp(
           dp_friction[n + 1] - height_ab/n/2*ambient.g*d[n],
           d[n],
           d_b,
           eta[n],
           eta_b,
           length/n/2,
           d_h,
           roughness,
           dp_small);
       else
         dp_friction[1] = WallFriction.pressureLoss_m_flow(
           m_flow[1],
           d_a,
           d[1],
           eta_a,
           eta[1],
           length/n/2,
           d_h,
           roughness,
           m_flow_small) + height_ab/n*ambient.g*d[1]/2;
         dp_friction[n] = WallFriction.pressureLoss_m_flow(
           m_flow[n + 1],
           d[n],
           d_b,
           eta[n],
           eta_b,
           length/n/2,
           d_h,
           roughness,
           m_flow_small) + height_ab/n/2*ambient.g*d[n];
       end if;
        
       annotation (
         Icon(
           Rectangle(extent=[-24,40; 24,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[-72,40; -24,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[24,40; 72,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[72,40; 96,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[-96,40; -72,-40], style(color=0, rgbcolor={0,0,0}))),
         Diagram,
         Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
 
</html>",      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
        
     end PartialDistributedFlow;
      
   end Symmetric;
    
   package Asymmetric 
      "Modelica_Fluid inspired discretization using the finite volume method on a staggered grid with static momentum balance (asymmetric)" 
     //extends Interfaces.PartialFluidDiscretization(isSymmetric=false);
      
      extends BaseClasses.FiniteVolumeOnStaggeredGridWithCorrelations(isSymmetric=
            false);
      
     redeclare replaceable partial model extends PartialDistributedFlow 
        
       // Taken one-to-one from Modelica_Fluid
        
       // Covered herein:
       //   TBD
        
     equation 
       //Static momentum balance
       0 = - (medium[1].p - port_a_p) - dp_friction[1];
       for i in 2:n loop
         0 = - (medium[i].p - medium[i - 1].p) - dp_friction[i];
       end for;
       0 = - (port_b_p - medium[n].p) - dp_friction[n + 1];
        
       // Source terms
       Qs_flow = heatTransfer.Q_flow;
       Ws_flow = zeros(n);
       ms_flow = zeros(n);
       msXi_flow = zeros(n, Medium.nXi);
       Vi = ones(n)*V/n;
       //Pressure drop and gravity
       if from_dp and not WallFriction.dp_is_zero then
         dp_friction[1] = 0 
            "Asymmetric model without momentum balance between port_a and the first control volume";
         m_flow[n + 1] = WallFriction.massFlowRate_dp(
           dp_friction[n + 1] - height_ab/n*ambient.g*d[n],
           d[n],
           d_b,
           eta[n],
           eta_b,
           length/n,
           d_h,
           roughness,
           dp_small);
       else
         dp_friction[1] = 0 
            "Asymmetric model without momentum balance between port_a and the first control volume";
         dp_friction[n] = WallFriction.pressureLoss_m_flow(
           m_flow[n + 1],
           d[n],
           d_b,
           eta[n],
           eta_b,
           length/n,
           d_h,
           roughness,
           m_flow_small) + height_ab/n*ambient.g*d[n];
       end if;
        
       annotation (
         Icon(
           Rectangle(extent=[-24,40; 24,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[-72,40; -24,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[24,40; 72,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[72,40; 120,-40], style(color=0, rgbcolor={0,0,0}))),
         Diagram,
         Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
 
</html>",      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
     end PartialDistributedFlow;
      
   end Asymmetric;
 end StaticMomentum;
  
 package DynamicMomentum "Discretizations with dynamic momentum balances" 
  extends Icons.VariantLibrary;
    
   package Asymmetric 
      "Discretization using the finite volume method on a staggered grid with dynamic momentum balance (asymmetric)" 
      
     extends BaseClasses.FiniteVolumeOnStaggeredGridWithCorrelations(isSymmetric=
           false);
      
     redeclare replaceable partial model extends PartialDistributedFlow 
        
       SI.MomentumFlux[n] I_flow 
          "Momentum flux of fluid through cv (staggered momentum balance)";
       SI.MomentumFlux I_flow_a "Momentum flux of fluid across port_a";
       SI.MomentumFlux I_flow_b "Momentum flux of fluid across port_b";
        
       // Covered herein:
       //   TBD
        
     equation 
       0 = port_a_p - medium[1].p 
          "Asymmetric model, no momentum balance to the 'left' of the first volume";
       // Dynamic momentum balance
       for i in 2:n loop
         der(m_flow[i]) * length/n = I_flow[i-1] - I_flow[i] - area * (medium[i].p - medium[i-1].p) - area * dp_friction[i];
       end for;
       der(m_flow[n+1])*length/n = I_flow[n] - I_flow_b - area * (port_b_p - medium[n].p) - area * dp_friction[n+1];
        
       // Central difference scheme for momentum flux calculation
       I_flow_a = 0 "Unused in this asymmetrical model";
       for i in 1:n loop
         I_flow[i]=(m_flow[i]+m_flow[i+1])^2/(4*medium[i].d*area);
       end for;
       I_flow_b=m_flow[n+1]^2/(d_b*area);
        
       // Source terms
       Qs_flow = heatTransfer.Q_flow;
       Ws_flow = zeros(n);
       ms_flow = zeros(n);
       msXi_flow = zeros(n, Medium.nXi);
       Vi = ones(n)*V/n;
       //Pressure drop and gravity
       if from_dp and not WallFriction.dp_is_zero then
         dp_friction[1] = 0 
            "Asymmetric model without momentum balance between port_a and the first control volume";
         m_flow[n + 1] = WallFriction.massFlowRate_dp(
           dp_friction[n + 1] - height_ab/n*ambient.g*d[n],
           d[n],
           d_b,
           eta[n],
           eta_b,
           length/n,
           d_h,
           roughness,
           dp_small);
       else
         dp_friction[1] = 0 
            "Asymmetric model without momentum balance between port_a and the first control volume";
         dp_friction[n] = WallFriction.pressureLoss_m_flow(
           m_flow[n + 1],
           d[n],
           d_b,
           eta[n],
           eta_b,
           length/n,
           d_h,
           roughness,
           m_flow_small) + height_ab/n*ambient.g*d[n];
       end if;
        
       annotation (
         Icon(
           Rectangle(extent=[-24,40; 24,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[-72,40; -24,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[24,40; 72,-40], style(color=0, rgbcolor={0,0,0})),
           Rectangle(extent=[72,40; 120,-40], style(color=0, rgbcolor={0,0,0}))),
         Diagram,
         Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
 
</html>",      revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
     end PartialDistributedFlow;
      
   end Asymmetric;
 end DynamicMomentum;
end FiniteVolumeOnStaggeredGrid;
