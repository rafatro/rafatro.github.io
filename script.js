const norm_s_inv_answers = {
  p75000:0.6744897502,
  p80000:0.8416212327,
  p85000:1.0364333907,
  p87500:1.1503493805,
  p90000:1.2815515641,
  p92500:1.4395314724,
  p95000:1.6448536251,
  p97500:1.9599639861,
  p98000:2.0537489090,
  p98750:2.2414027265,
  p99000:2.3263478744,
  p99500:2.5758293064,
  p99750:2.8070337701,
  p99900:3.0902323047,
  p99950:3.2905267283,
  p99990:3.7190164821,
  p99995:3.8905918843
};

function norm_s_inv(prob){
  return(norm_s_inv_answers["p"+prob*100000])
};

function calculateMainPhrase(){
  //Hide everything:
  document.getElementById("ssp").style.display="none";
  document.getElementById("sspf").style.display="none";
  document.getElementById("ssv").style.display="none";
  document.getElementById("ssvf").style.display="none";
  document.getElementById("2pt").style.display="none";
  document.getElementById("2pf").style.display="none";
  document.getElementById("2vt").style.display="none";
  document.getElementById("2vf").style.display="none";
  document.getElementById("FAQ_title").style.display="none";
  document.getElementById("FAQ_alpha").style.display="none";
  document.getElementById("FAQ_rate").style.display="none";
  document.getElementById("FAQ_defineInterval").style.display="none";
  document.getElementById("FAQ_rateForFixed").style.display="none";
  document.getElementById("FAQ_calc_ssp").style.display="none";
  document.getElementById("FAQ_calc_sspf").style.display="none";
  document.getElementById("FAQ_calc_ssv").style.display="none";
  document.getElementById("FAQ_calc_ssvf").style.display="none";

  if (document.getElementById("samples").value=="1") {
    document.getElementById("textOp1").innerText="";
    document.getElementById("textOp2").innerText="";
    document.getElementById("optionCI").style.display="block";
    document.getElementById("optionStatsig").style.display="none";
    if (document.getElementById("goal").value=="statsig") document.getElementById("goal").value="confint";
  } else {
    document.getElementById("textOp1").innerText="of each one";
    document.getElementById("textOp2").innerText="for the difference between them";
    document.getElementById("optionStatsig").style.display="block";
    document.getElementById("optionCI").style.display="none";
    if (document.getElementById("goal").value=="confint") document.getElementById("goal").value="statsig";
  }
  var testType="";
  if (document.getElementById("samples").value=="1"){
    if (document.getElementById("metric").value=="rate") document.getElementById("textOp3").innerText="different from a given rate";
    else document.getElementById("textOp3").innerText="different from a given value";
  }
  else {
    if (document.getElementById("metric").value=="rate") document.getElementById("textOp3").innerText="larger than a given p.p.";
    else document.getElementById("textOp3").innerText="larger than a given value";
  }
}

function doMainPhrase(){
  if (document.getElementById("samples").value=="1"){
    if (document.getElementById("metric").value=="rate"){
      if (document.getElementById("goal").value=="confint"){
        document.getElementById("ssp").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        //document.getElementById("FAQ_rate").style.display="block";
        //document.getElementById("FAQ_defineInterval").style.display="block";
        //document.getElementById("FAQ_calc_ssp").style.display="block";
        calculateOneSampleRate();
      } else {
        document.getElementById("sspf").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        //document.getElementById("FAQ_rateForFixed").style.display="block";
        //document.getElementById("FAQ_calc_sspf").style.display="block";
        calculateOneSampleRateFixed();
      }
    }
    else {
      if (document.getElementById("goal").value=="confint") {
        document.getElementById("ssv").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        //document.getElementById("FAQ_calc_ssv").style.display="block";
        //add more
        calculateOneSampleAvg();
      } else {
        document.getElementById("ssvf").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        //document.getElementById("FAQ_calc_ssvf").style.display="block";
        // add more
        calculateOneSampleAvgFixed();
      }
    }
  }
  else {
    if (document.getElementById("metric").value=="rate"){
      if (document.getElementById("goal").value=="statsig") {
        document.getElementById("2pt").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        // add more
        calculate2SampleRate();
      }
      else {
        document.getElementById("2pf").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        // add more
        calculate2SampleRateF();
      }
    }
    else {
      if (document.getElementById("goal").value=="statsig") {
        document.getElementById("2vt").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        // add more
        calculate2SampleAvg();
      }
      else {
        document.getElementById("2vf").style.display="block";
        //document.getElementById("FAQ_title").style.display="block";
        //document.getElementById("FAQ_alpha").style.display="block";
        // add more
        calculate2SampleAvgF();
      }
    }
  }
}

function adjustInputFieldSize(inputFieldObj){
  inputFieldObj.style.width=Math.max(1,inputFieldObj.value.length)+2+"ch";
  console.log(inputFieldObj.style.width)
}

function calculateOneSampleRate(){
  //document.getElementById("oneSampleRateP").size = Math.max(1,document.getElementById("oneSampleRateP").value.length);
  adjustInputFieldSize(document.getElementById("oneSampleRateP"))
  document.getElementById("oneSampleRateCI").size = Math.max(1,document.getElementById("oneSampleRateCI").value.length);
  // add validation oneSampleRateRateP !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleRateP").value.split(".").length == 1){
    var p_precision=0;
  } else {
    var p_precision=document.getElementById("oneSampleRateP").value.split(".")[1].length;
  }
  var p = parseFloat(document.getElementById("oneSampleRateP").value);
  // add validation!!! tem que ser número between 0 and 100%
  // add validation oneSampleRateCI !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleRateCI").value.split(".").length == 1){
    var ci_precision=0;
  } else {
    var ci_precision=document.getElementById("oneSampleRateCI").value.split(".")[1].length;
  }
  var precision = Math.max(p_precision,ci_precision);
  var ci = parseFloat(document.getElementById("oneSampleRateCI").value);
  // add validation!!! tem que ser número positive non zero, under 100 p.p.
  p = p/100;
  ci = ci/100;
  document.getElementById("oneSampleRateMin").innerText=((p-ci)*100).toFixed(precision)+"%";
  document.getElementById("oneSampleRateMax").innerText=((p+ci)*100).toFixed(precision)+"%";
  document.getElementById("oneSampleRateHyp1").innerText=(p*100).toFixed(precision)+"%";
  document.getElementById("oneSampleRateHyp2").innerText=(p*100).toFixed(precision)+"%";
  document.getElementById("oneSampleRateHyp3").innerText=(p*100).toFixed(precision)+"%";
  document.getElementById("oneSampleRateHyp4").innerText=(p*100).toFixed(precision)+"%";
  var stdev = Math.sqrt(p*(1-p));
  var conf = parseFloat(document.getElementById("oneSampleRateConf").value);
  var z = norm_s_inv(conf+(1-conf)/2); //this test is two-tailed
  var n = Math.pow(z*stdev/ci,2);
  var nPower = Math.pow((z+norm_s_inv(0.8))*stdev/ci,2);
  document.getElementById("oneSampleRateN").innerText=Math.ceil(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  document.getElementById("oneSampleRateNPower").innerText=Math.ceil(nPower).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  draw_oneSampleRateCI(conf,ci,p,precision);
};

function calculateOneSampleRateFixed(){
  // Clear error mesages:
  document.getElementById("oneSampleRateFixedError1").style.display='none';
  document.getElementById("oneSampleRateFixedError2").style.display='none';
  //
  document.getElementById("oneSampleRateFixedP").size = Math.max(1,document.getElementById("oneSampleRateFixedP").value.length);
  document.getElementById("oneSampleRateFixedSampled").size = Math.max(1,document.getElementById("oneSampleRateFixedSampled").value.length);
  // add validations !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleRateFixedP").value.split(".").length == 1){
    var p_precision=0;
  } else {
    var p_precision=document.getElementById("oneSampleRateFixedP").value.split(".")[1].length;
  }
  var p = parseFloat(document.getElementById("oneSampleRateFixedP").value);
  // add validation!!! tem que ser número between 0 and 100%
  // add validation oneSampleRateFixedSampled !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleRateFixedSampled").value.split(".").length == 1){
    var s_precision=0;
  } else {
    var s_precision=document.getElementById("oneSampleRateFixedSampled").value.split(".")[1].length;
  }
  var s = parseFloat(document.getElementById("oneSampleRateFixedSampled").value);
  var precision = Math.max(p_precision,s_precision);
  // add validation!!! tem que ser número positive non zero, under 100 p.p.
  p = p/100;
  s = s/100;
  var delta = s-p;
  var stdev = Math.sqrt(p*(1-p));
  var conf = parseFloat(document.getElementById("oneSampleRateFixedConf").value);
  var type = document.getElementById("oneSampleRateFixedAns").value;
  document.getElementById("oneSampleRateFixedHyp1").innerText=(p*100).toFixed(precision)+"%";
  document.getElementById("oneSampleRateFixedHyp3").innerText=(p*100).toFixed(precision)+"%";
  document.getElementById("oneSampleRateFixedHyp4").innerText=(p*100).toFixed(precision)+"%";
  if (type=='not'){
    document.getElementById("oneSampleRateFixedTwoTailed").style.display='block';
    document.getElementById("oneSampleRateFixedRightTailed").style.display='none';
    document.getElementById("oneSampleRateFixedLeftTailed").style.display='none';
    document.getElementById("oneSampleRateFixedHyp2").innerText="different from "+(p*100).toFixed(precision)+"%";
    var z = norm_s_inv(conf+(1-conf)/2); //this test is two-tailed
    var n = Math.pow(z*stdev/Math.abs(delta),2);
    var nPower = Math.pow((z+norm_s_inv(0.8))*stdev/Math.abs(delta),2);
  } else {
    if (type=='not smaller') {
      document.getElementById("oneSampleRateFixedTwoTailed").style.display='none';
      document.getElementById("oneSampleRateFixedRightTailed").style.display='block';
      document.getElementById("oneSampleRateFixedLeftTailed").style.display='none';
      document.getElementById("oneSampleRateFixedHyp2").innerText="larger than "+(p*100).toFixed(precision)+"%";
    } else {
      document.getElementById("oneSampleRateFixedTwoTailed").style.display='none';
      document.getElementById("oneSampleRateFixedRightTailed").style.display='none';
      document.getElementById("oneSampleRateFixedLeftTailed").style.display='block';
      document.getElementById("oneSampleRateFixedHyp2").innerText="smaller than "+(p*100).toFixed(precision)+"%";
    }
    if (type=='not smaller' && delta<0){
      n, nPower = Infinity
      document.getElementById("oneSampleRateFixedError1").style.display='block';
    } else if (type=='not larger' && delta>0){
      n, nPower = Infinity
      document.getElementById("oneSampleRateFixedError2").style.display='block';
    } else {
      var z = norm_s_inv(conf); //this test is one-tailed
      var n = Math.pow(z*stdev/delta,2);
      var nPower = Math.pow((z+norm_s_inv(0.8))*stdev/delta,2);
    }
  }
  document.getElementById("oneSampleRateFixedN").innerText=Math.ceil(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  document.getElementById("oneSampleRateFixedNPower").innerText=Math.ceil(nPower).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  draw_oneSampleRateFixedS(conf,s,p,precision,type);
};

function calculateOneSampleAvg(){
  document.getElementById("oneSampleAvgStdev").size = Math.max(1,document.getElementById("oneSampleAvgStdev").value.length);
  document.getElementById("oneSampleAvgCI").size = Math.max(1,document.getElementById("oneSampleAvgCI").value.length);
  // add validation !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleAvgStdev").value.split(".").length == 1){
    var stdev_precision=0;
  } else {
    var stdev_precision=document.getElementById("oneSampleAvgStdev").value.split(".")[1].length;
  }
  // add validation!!!
  // add validation oneSampleAvgCI !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleAvgCI").value.split(".").length == 1){
    var ci_precision=0;
  } else {
    var ci_precision=document.getElementById("oneSampleAvgCI").value.split(".")[1].length;
  }
  var precision = Math.max(stdev_precision,ci_precision);
  var ci = parseFloat(document.getElementById("oneSampleAvgCI").value);
  // add validation!!! tem que ser número positive non zero
  var stdev = parseFloat(document.getElementById("oneSampleAvgStdev").value);
  var conf = parseFloat(document.getElementById("oneSampleAvgConf").value);
  var z = norm_s_inv(conf+(1-conf)/2); //this test is two-tailed
  var n = Math.pow(z*stdev/ci,2);
  document.getElementById("oneSampleAvgN").innerText=Math.ceil(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  //draw_oneSampleRateCI(conf,ci,p,precision);
};

function calculateOneSampleAvgFixed(){
  // Clear error mesages:
  document.getElementById("oneSampleAvgFixedError1").style.display='none';
  document.getElementById("oneSampleAvgFixedError2").style.display='none';
  //
  document.getElementById("oneSampleAvgFixedStdev").size = Math.max(1,document.getElementById("oneSampleAvgFixedStdev").value.length);
  document.getElementById("oneSampleAvgFixedAvg").size = Math.max(1,document.getElementById("oneSampleAvgFixedAvg").value.length);
  document.getElementById("oneSampleAvgFixedSampled").size = Math.max(1,document.getElementById("oneSampleAvgFixedSampled").value.length);
  // add validations !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleAvgFixedAvg").value.split(".").length == 1){
    var avg_precision=0;
  } else {
    var avg_precision=document.getElementById("oneSampleAvgFixedAvg").value.split(".")[1].length;
  }
  var avg = parseFloat(document.getElementById("oneSampleAvgFixedAvg").value);
  // add validation!!! tem que ser número between 0 and 100%
  // add validation oneSampleRateFixedSampled !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("oneSampleAvgFixedSampled").value.split(".").length == 1){
    var s_precision=0;
  } else {
    var s_precision=document.getElementById("oneSampleAvgFixedSampled").value.split(".")[1].length;
  }
  var s = parseFloat(document.getElementById("oneSampleAvgFixedSampled").value);
  var precision = Math.max(avg_precision,s_precision);
  // add validation!!! tem que ser número positive non zero, under 100 p.p.
  var delta = s-avg;
  var stdev = parseFloat(document.getElementById("oneSampleAvgFixedStdev").value);
  var conf = parseFloat(document.getElementById("oneSampleAvgFixedConf").value);
  document.getElementById("oneSampleAvgFixedHyp1").innerText=avg.toFixed(precision);
  document.getElementById("oneSampleAvgFixedHyp3").innerText=avg.toFixed(precision);
  document.getElementById("oneSampleAvgFixedHyp4").innerText=avg.toFixed(precision);
  var type = document.getElementById("oneSampleAvgFixedAns").value;
  if (type=='not'){
    document.getElementById("oneSampleAvgFixedTwoTailed").style.display='block';
    document.getElementById("oneSampleAvgFixedRightTailed").style.display='none';
    document.getElementById("oneSampleAvgFixedLeftTailed").style.display='none';
    document.getElementById("oneSampleAvgFixedHyp2").innerText="different from "+avg.toFixed(precision);
    var z = norm_s_inv(conf+(1-conf)/2); //this test is two-tailed
    var n = Math.pow(z*stdev/Math.abs(delta),2);
    var nPower = Math.pow((z+norm_s_inv(0.8))*stdev/Math.abs(delta),2);
  } else {
    if (type=='not smaller') {
      document.getElementById("oneSampleAvgFixedTwoTailed").style.display='none';
      document.getElementById("oneSampleAvgFixedRightTailed").style.display='block';
      document.getElementById("oneSampleAvgFixedLeftTailed").style.display='none';
      document.getElementById("oneSampleAvgFixedHyp2").innerText="larger than "+avg.toFixed(precision);
    } else {
      document.getElementById("oneSampleAvgFixedTwoTailed").style.display='none';
      document.getElementById("oneSampleAvgFixedRightTailed").style.display='none';
      document.getElementById("oneSampleAvgFixedLeftTailed").style.display='block';
      document.getElementById("oneSampleAvgFixedHyp2").innerText="smaller than "+avg.toFixed(precision);
    }
    if (type=='not smaller' && delta<0){
      n, nPower = Infinity
      document.getElementById("oneSampleAvgFixedError1").style.display='block';
    } else if (type=='not larger' && delta>0){
      n, nPower = Infinity
      document.getElementById("oneSampleAvgFixedError2").style.display='block';
    } else {
      var z = norm_s_inv(conf); //this test is one-tailed
      var n = Math.pow(z*stdev/delta,2);
      var nPower = Math.pow((z+norm_s_inv(0.8))*stdev/delta,2);
    }
  }
  document.getElementById("oneSampleAvgFixedN").innerText=Math.ceil(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  document.getElementById("oneSampleAvgFixedNPower").innerText=Math.ceil(nPower).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  //draw_oneSampleAvgFixedS(conf,s,p,precision,type);
};

function calculate2SampleRate(){
  // Clear error mesages:
  document.getElementById("2SampleRateError1").style.display='none';
  document.getElementById("2SampleRateError2").style.display='none';

  document.getElementById("2SampleRateP1").size = Math.max(1,document.getElementById("2SampleRateP1").value.length);
  document.getElementById("2SampleRateP2").size = Math.max(1,document.getElementById("2SampleRateP2").value.length);
  // add validation oneSampleRateRateP !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("2SampleRateP1").value.split(".").length == 1){
    var p1_precision=0;
  } else {
    var p1_precision=document.getElementById("2SampleRateP1").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleRateP2").value.split(".").length == 1){
    var p2_precision=0;
  } else {
    var p2_precision=document.getElementById("2SampleRateP2").value.split(".")[1].length;
  }
  var p1 = parseFloat(document.getElementById("2SampleRateP1").value);
  var p2 = parseFloat(document.getElementById("2SampleRateP2").value);
  // add validation!!! tem que ser número between 0 and 100%
  // add validation oneSampleRateCI !!! tem que ter sõ número e até 1 ponto
  var precision = Math.max(p1_precision,p2_precision);
  p1 = p1/100;
  p2 = p2/100;
  delta = p2-p1;
  document.getElementById("2SampleRateDelta").innerText=(delta*100).toFixed(precision)+" p.p.";
  var n_proportion = parseFloat(document.getElementById("2SampleRateSplit").value)/100; // proportion of n1 inside n1+n2
  var p_pooled = n_proportion*p1+(1-n_proportion)*p2;
  var conf = parseFloat(document.getElementById("2SampleRateConf").value);
  var type = document.getElementById("2SampleRateAns").value;
  if (type == "not equal") {
    document.getElementById("2SampleRateTwoTailed").style.display='block';
    document.getElementById("2SampleRateRightTailed").style.display='none';
    document.getElementById("2SampleRateLeftTailed").style.display='none';
    document.getElementById("2SampleRateHyp1").innerText="Samples come from populations with different rates";
    var z = norm_s_inv(conf+(1-conf)/2); //this test is two-tailed
    var n1 = Math.pow(z*Math.sqrt(p_pooled*(1-p_pooled))/delta,2)/(1-n_proportion);
    var n2 = n1*(1/n_proportion-1);
    var n1Power = Math.pow((z+norm_s_inv(0.8))*Math.sqrt(p_pooled*(1-p_pooled))/delta,2)/(1-n_proportion);
    var n2Power = n1Power*(1/n_proportion-1);
  } else {
    if (type=='not smaller') {
      document.getElementById("2SampleRateTwoTailed").style.display='none';
      document.getElementById("2SampleRateRightTailed").style.display='block';
      document.getElementById("2SampleRateLeftTailed").style.display='none';
      document.getElementById("2SampleRateHyp1").innerText="The second sample comes from a population with a larger rate than the first";
    } else {
      document.getElementById("2SampleRateTwoTailed").style.display='none';
      document.getElementById("2SampleRateRightTailed").style.display='none';
      document.getElementById("2SampleRateLeftTailed").style.display='block';
      document.getElementById("2SampleRateHyp1").innerText="The second sample comes from a population with a smaller rate than the first";
    }
    if (type=='not smaller' && delta<0){
      var n1, n2, n1Power, n2Power = Infinity
      document.getElementById("2SampleRateError1").style.display='block';
    } else if (type=='not larger' && delta>0){
      var n1, n2 = Infinity
      var n1, n2, n1Power, n2Power = Infinity
      document.getElementById("2SampleRateError2").style.display='block';
    } else {
      var z = norm_s_inv(conf); //this test is one-tailed
      var n1 = Math.pow(z*Math.sqrt(p_pooled*(1-p_pooled))/delta,2)/(1-n_proportion);
      var n2 = n1*(1/n_proportion-1);
      var n1Power = Math.pow((z+norm_s_inv(0.8))*Math.sqrt(p_pooled*(1-p_pooled))/delta,2)/(1-n_proportion);
      var n2Power = n1Power*(1/n_proportion-1);
    }
  } // I have to test more against online calculators... I just ound 1 that performs unequal sample size and one-tail (http://www2.psych.purdue.edu/~gfrancis/calculators/proportion_test_two_sample.shtml)
  if (n_proportion == 0.5){
    document.getElementById("2SampleRateN_phrase").innerText="2 samples of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
    document.getElementById("2SampleRateNPower_phrase").innerText="2 samples of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
  }
  else {
    document.getElementById("2SampleRateN_phrase").innerText="Sample sizes of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
    document.getElementById("2SampleRateNPower_phrase").innerText="Sample sizes of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
  }
  //draw_oneSampleRateCI(conf,ci,p,precision);
};

function calculate2SampleRateF(){
  // Clear error mesages:
  document.getElementById("2SampleRateFError1").style.display='none';
  document.getElementById("2SampleRateFError2").style.display='none';
  document.getElementById("2SampleRateFError3").style.display='none';
  document.getElementById("2SampleRateFError4").style.display='none';

  document.getElementById("2SampleRateFP1").size = Math.max(1,document.getElementById("2SampleRateFP1").value.length);
  document.getElementById("2SampleRateFP2").size = Math.max(1,document.getElementById("2SampleRateFP2").value.length);
  document.getElementById("2SampleRateFMinDelta").size = Math.max(1,document.getElementById("2SampleRateFMinDelta").value.length);
  // add validation oneSampleRateRateP !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("2SampleRateFP1").value.split(".").length == 1){
    var p1_precision=0;
  } else {
    var p1_precision=document.getElementById("2SampleRateFP1").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleRateFP2").value.split(".").length == 1){
    var p2_precision=0;
  } else {
    var p2_precision=document.getElementById("2SampleRateFP2").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleRateFMinDelta").value.split(".").length == 1){
    var ci_precision=0;
  } else {
    var ci_precision=document.getElementById("2SampleRateFMinDelta").value.split(".")[1].length;
  }
  var p1 = parseFloat(document.getElementById("2SampleRateFP1").value);
  var p2 = parseFloat(document.getElementById("2SampleRateFP2").value);
  // add validation!!! tem que ser número between 0 and 100%
  // add validation oneSampleRateCI !!! tem que ter sõ número e até 1 ponto
  var precision = Math.max(p1_precision,p2_precision,ci_precision);
  p1 = p1/100;
  p2 = p2/100;
  delta = p2-p1;
  document.getElementById("2SampleRateFDelta").innerText=(delta*100).toFixed(precision)+" p.p.";
  document.getElementById("2SampleRateFDelta_error2").innerText=(delta*100).toFixed(precision)+" p.p.";
  document.getElementById("2SampleRateFDelta_error4").innerText=(delta*100).toFixed(precision)+" p.p.";
  var minDelta = parseFloat(document.getElementById("2SampleRateFMinDelta").value);
  minDelta = minDelta/100;
  document.getElementById("2SampleRateFMinDelta_error2").innerText=(minDelta*100).toFixed(precision)+" p.p.";
  document.getElementById("2SampleRateFMinDelta_error4").innerText=(minDelta*100).toFixed(precision)+" p.p.";
  document.getElementById("2SampleRateFHyp1").innerText=(minDelta*100).toFixed(precision);
  document.getElementById("2SampleRateFHyp2").innerText=(minDelta*100).toFixed(precision);
  document.getElementById("2SampleRateFHyp3").innerText=(minDelta*100).toFixed(precision);
  document.getElementById("2SampleRateFHyp4").innerText=(minDelta*100).toFixed(precision);
  var n_proportion = parseFloat(document.getElementById("2SampleRateFSplit").value)/100; // proportion of n1 inside n1+n2
  var conf = parseFloat(document.getElementById("2SampleRateFConf").value);
  if (delta>0 && minDelta<0) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleRateFError1").style.display='block';
  } else if (delta>0 && minDelta>=delta) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleRateFError2").style.display='block';
  } else if (delta<0 && minDelta>0) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleRateFError3").style.display='block';
  } else if (delta<0 && minDelta<=delta) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleRateFError4").style.display='block';
  }
  else {
    var z = norm_s_inv(conf); //this test is one-tailed
    var n1 = Math.pow(z*Math.sqrt(p1*(1-p1)+p2*(1-p2)/(1/n_proportion-1))/(Math.abs(delta)-Math.abs(minDelta)),2);
    var n2 = n1*(1/n_proportion-1);
    var n1Power = Math.pow((z+norm_s_inv(0.8))*Math.sqrt(p1*(1-p1)+p2*(1-p2)/(1/n_proportion-1))/(Math.abs(delta)-Math.abs(minDelta)),2);
    var n2Power = n1Power*(1/n_proportion-1);
  }
  // I have to test more against online calculators... I didn~t find any
  if (n_proportion == 0.5) {
    document.getElementById("2SampleRateFN_phrase").innerText="2 samples of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
    document.getElementById("2SampleRateFNPower_phrase").innerText="2 samples of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
  }
  else {
    document.getElementById("2SampleRateFN_phrase").innerText="Sample sizes of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
    document.getElementById("2SampleRateFNPower_phrase").innerText="Sample sizes of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
  }
  //draw_oneSampleRateCI(conf,ci,p,precision);
};

function calculate2SampleAvg(){
  // Clear error mesages:
  document.getElementById("2SampleAvgError1").style.display='none';
  document.getElementById("2SampleAvgError2").style.display='none';

  document.getElementById("2SampleAvg1").size = Math.max(1,document.getElementById("2SampleAvg1").value.length);
  document.getElementById("2SampleAvgStdev1").size = Math.max(1,document.getElementById("2SampleAvgStdev1").value.length);
  document.getElementById("2SampleAvg2").size = Math.max(1,document.getElementById("2SampleAvg2").value.length);
  document.getElementById("2SampleAvgStdev2").size = Math.max(1,document.getElementById("2SampleAvgStdev2").value.length);
  // add validation oneSampleRateRateP !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("2SampleAvg1").value.split(".").length == 1){
    var avg1_precision=0;
  } else {
    var avg1_precision=document.getElementById("2SampleAvg1").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleAvg2").value.split(".").length == 1){
    var avg2_precision=0;
  } else {
    var avg2_precision=document.getElementById("2SampleAvg2").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleAvgStdev1").value.split(".").length == 1){
    var stdev1_precision=0;
  } else {
    var stdev1_precision=document.getElementById("2SampleAvgStdev1").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleAvgStdev2").value.split(".").length == 1){
    var stdev2_precision=0;
  } else {
    var stdev2_precision=document.getElementById("2SampleAvgStdev2").value.split(".")[1].length;
  }
  var avg1 = parseFloat(document.getElementById("2SampleAvg1").value);
  var avg2 = parseFloat(document.getElementById("2SampleAvg2").value);
  var stdev1 = parseFloat(document.getElementById("2SampleAvgStdev1").value);
  var stdev2 = parseFloat(document.getElementById("2SampleAvgStdev2").value);
  // add validation!!! tem que ser número between 0 and 100%
  // add validation oneSampleRateCI !!! tem que ter sõ número e até 1 ponto
  var precision = Math.max(avg1_precision,avg2_precision,stdev1_precision,stdev2_precision);
  delta = avg2-avg1;
  document.getElementById("2SampleAvgDelta").innerText=(delta).toFixed(precision);
  var n_proportion = parseFloat(document.getElementById("2SampleAvgSplit").value)/100; // proportion of n1 inside n1+n2
  var conf = parseFloat(document.getElementById("2SampleAvgConf").value);
  var type = document.getElementById("2SampleAvgAns").value;
  if (type == "not equal") {
    document.getElementById("2SampleAvgTwoTailed").style.display='block';
    document.getElementById("2SampleAvgRightTailed").style.display='none';
    document.getElementById("2SampleAvgLeftTailed").style.display='none';
    document.getElementById("2SampleAvgHyp1").innerText="Samples come from populations with different mean values";
    var z = norm_s_inv(conf+(1-conf)/2); //this test is two-tailed
    var n1 = Math.pow(z*Math.sqrt(Math.pow(stdev1,2)+Math.pow(stdev2,2)/(1/n_proportion-1))/Math.abs(delta),2);
    var n2 = n1*(1/n_proportion-1);
    var n1Power = Math.pow((z+norm_s_inv(0.8))*Math.sqrt(Math.pow(stdev1,2)+Math.pow(stdev2,2)/(1/n_proportion-1))/Math.abs(delta),2);
    var n2Power = n1Power*(1/n_proportion-1);
  } else {
    if (type=='not smaller') {
      document.getElementById("2SampleAvgTwoTailed").style.display='none';
      document.getElementById("2SampleAvgRightTailed").style.display='block';
      document.getElementById("2SampleAvgLeftTailed").style.display='none';
      document.getElementById("2SampleAvgHyp1").innerText="The second sample comes from a population with a larger mean value than the first";
    } else {
      document.getElementById("2SampleAvgTwoTailed").style.display='none';
      document.getElementById("2SampleAvgRightTailed").style.display='none';
      document.getElementById("2SampleAvgLeftTailed").style.display='block';
      document.getElementById("2SampleAvgHyp1").innerText="The second sample comes from a population with a smaller mean value than the first";
    }
    if (type=='not smaller' && delta<0){
      var n1, n2, n1Power, n2Power = Infinity
      document.getElementById("2SampleAvgError1").style.display='block';
    } else if (type=='not larger' && delta>0){
      var n1, n2, n1Power, n2Power = Infinity
      document.getElementById("2SampleAvgError2").style.display='block';
    } else {
      var z = norm_s_inv(conf); //this test is one-tailed
      var n1 = Math.pow(z*Math.sqrt(Math.pow(stdev1,2)+Math.pow(stdev2,2)/(1/n_proportion-1))/Math.abs(delta),2);
      var n2 = n1*(1/n_proportion-1);
      var n1Power = Math.pow((z+norm_s_inv(0.8))*Math.sqrt(Math.pow(stdev1,2)+Math.pow(stdev2,2)/(1/n_proportion-1))/Math.abs(delta),2);
      var n2Power = n1Power*(1/n_proportion-1);
    }
  }
  if (n_proportion == 0.5) {
    document.getElementById("2SampleAvgN_phrase").innerText="2 samples of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
    document.getElementById("2SampleAvgNPower_phrase").innerText="2 samples of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
  }
  else {
    document.getElementById("2SampleAvgN_phrase").innerText="Sample sizes of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
    document.getElementById("2SampleAvgNPower_phrase").innerText="Sample sizes of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
  }
  //draw_oneSampleRateCI(conf,ci,p,precision);
};

function calculate2SampleAvgF(){
  // Clear error mesages:
  document.getElementById("2SampleAvgFError1").style.display='none';
  document.getElementById("2SampleAvgFError2").style.display='none';
  document.getElementById("2SampleAvgFError3").style.display='none';
  document.getElementById("2SampleAvgFError4").style.display='none';

  document.getElementById("2SampleAvgF1").size = Math.max(1,document.getElementById("2SampleAvgF1").value.length);
  document.getElementById("2SampleAvgFStdev1").size = Math.max(1,document.getElementById("2SampleAvgFStdev1").value.length);
  document.getElementById("2SampleAvgF2").size = Math.max(1,document.getElementById("2SampleAvgF2").value.length);
  document.getElementById("2SampleAvgFStdev2").size = Math.max(1,document.getElementById("2SampleAvgFStdev2").value.length);
  document.getElementById("2SampleAvgFMinDelta").size = Math.max(1,document.getElementById("2SampleAvgFMinDelta").value.length);
  // add validation oneSampleRateRateP !!! tem que ter sõ número e até 1 ponto
  if (document.getElementById("2SampleAvgF1").value.split(".").length == 1){
    var avg1_precision=0;
  } else {
    var avg1_precision=document.getElementById("2SampleAvgF1").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleAvgF2").value.split(".").length == 1){
    var avg2_precision=0;
  } else {
    var avg2_precision=document.getElementById("2SampleAvgF2").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleAvgFMinDelta").value.split(".").length == 1){
    var minDelta_precision=0;
  } else {
    var minDelta_precision=document.getElementById("2SampleAvgFMinDelta").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleAvgFStdev1").value.split(".").length == 1){
    var stdev1_precision=0;
  } else {
    var stdev1_precision=document.getElementById("2SampleAvgFStdev1").value.split(".")[1].length;
  }
  if (document.getElementById("2SampleAvgFStdev2").value.split(".").length == 1){
    var stdev2_precision=0;
  } else {
    var stdev2_precision=document.getElementById("2SampleAvgFStdev2").value.split(".")[1].length;
  }
  var avg1 = parseFloat(document.getElementById("2SampleAvgF1").value);
  var avg2 = parseFloat(document.getElementById("2SampleAvgF2").value);
  var stdev1 = parseFloat(document.getElementById("2SampleAvgFStdev1").value);
  var stdev2 = parseFloat(document.getElementById("2SampleAvgFStdev2").value);
  // add validation!!! tem que ser número between 0 and 100%
  // add validation oneSampleRateCI !!! tem que ter sõ número e até 1 ponto
  var precision = Math.max(avg1_precision,avg2_precision,minDelta_precision,stdev1_precision,stdev2_precision);
  delta = avg2 - avg1;
  document.getElementById("2SampleAvgFDelta").innerText=delta.toFixed(precision);
  document.getElementById("2SampleAvgFDelta_error2").innerText=delta.toFixed(precision);
  document.getElementById("2SampleAvgFDelta_error4").innerText=delta.toFixed(precision);
  var minDelta = parseFloat(document.getElementById("2SampleAvgFMinDelta").value);
  document.getElementById("2SampleAvgFMinDelta_error2").innerText=minDelta.toFixed(precision);
  document.getElementById("2SampleAvgFMinDelta_error4").innerText=minDelta.toFixed(precision);
  document.getElementById("2SampleAvgFHyp1").innerText=minDelta.toFixed(precision);
  document.getElementById("2SampleAvgFHyp2").innerText=minDelta.toFixed(precision);
  document.getElementById("2SampleAvgFHyp3").innerText=minDelta.toFixed(precision);
  document.getElementById("2SampleAvgFHyp4").innerText=minDelta.toFixed(precision);
  var n_proportion = parseFloat(document.getElementById("2SampleAvgFSplit").value)/100; // proportion of n1 inside n1+n2
  var conf = parseFloat(document.getElementById("2SampleAvgFConf").value);
  if (delta>0 && minDelta<0) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleAvgFError1").style.display='block';
  } else if (delta>0 && minDelta>=delta) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleAvgFError2").style.display='block';
  } else if (delta<0 && minDelta>0) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleAvgFError3").style.display='block';
  } else if (delta<0 && minDelta<=delta) {
    var n1, n2, n1Power, n2Power = Infinity
    document.getElementById("2SampleAvgFError4").style.display='block';
  }
  else {
    var z = norm_s_inv(conf); //this test is one-tailed
    var n1 = Math.pow(z*Math.sqrt(Math.pow(stdev1,2)+Math.pow(stdev2,2)/(1/n_proportion-1))/(Math.abs(delta)-Math.abs(minDelta)),2);
    var n2 = n1*(1/n_proportion-1);
    var n1Power = Math.pow((z+norm_s_inv(0.8))*Math.sqrt(Math.pow(stdev1,2)+Math.pow(stdev2,2)/(1/n_proportion-1))/(Math.abs(delta)-Math.abs(minDelta)),2);
    var n2Power = n1Power*(1/n_proportion-1);
  }
  if (n_proportion == 0.5) {
    document.getElementById("2SampleAvgFN_phrase").innerText="2 samples of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
    document.getElementById("2SampleAvgFNPower_phrase").innerText="2 samples of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations each";
  }
  else {
    document.getElementById("2SampleAvgFN_phrase").innerText="Sample sizes of " + Math.ceil(n1).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
    document.getElementById("2SampleAvgFNPower_phrase").innerText="Sample sizes of " + Math.ceil(n1Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " and " + Math.ceil(n2Power).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") + " observations respectively";
  }
  //draw_oneSampleAvgFCI(conf,ci,p,precision);
};

function afterBodyLoad(){
  document.querySelectorAll('[data-bs-toggle="tooltip"]')
  .forEach(tooltip => {
    new bootstrap.Tooltip(tooltip)
  })
  calculateMainPhrase();
}