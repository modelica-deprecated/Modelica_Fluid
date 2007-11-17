model TestWaterPump2 
  import PC = Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics;
  function pumpConsChar2 = PC.quadraticPower (
    q_nom={0,0.001,0.0015}, W_nom={550,650,800});
  extends TestWaterPumpDefault(Pump1(M=0.1,redeclare function 
        powerCharacteristic = 
          pumpConsChar2,
      usePowerCharacteristic=true),
    Source(T=ambient.default_T_ambient),
    Sink(T=ambient.default_T_ambient));
  
  annotation (Diagram,
    experiment(StopTime=10),
    experimentSetupOutput);
  
end TestWaterPump2;
