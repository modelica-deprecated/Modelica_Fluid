model TestWaterPump2 
  import PC = Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics;
  function pumpConsChar2 = PC.quadraticPower (
    q_nom={0,0.001,0.0015}, W_nom={550,650,800});
  extends TestWaterPumpDefault(Pump1(M=0.1,redeclare function 
        powerCharacteristic = 
          pumpConsChar2), 
    Source(T=ambient.default_T_ambient), 
    SinkP1(T=ambient.default_T_ambient));
  
  annotation (Diagram);
  
end TestWaterPump2;