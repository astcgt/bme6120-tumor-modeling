/** Results sequence follows a hypothesis and its evidence, not file order. */
export function tellResults({slide,text,line,foot,table,math,scientificChart,curve,side,results}){
 const r=results.therapy;
 const percentiles=(name,c,color,key='total',transform=x=>x)=>[
  curve((name+' P10').trim(),c[key+'_p10'].map(transform),color,true),
  curve((name+' 中位').trim(),c[key+'_median'].map(transform),color),
  curve((name+' P90').trim(),c[key+'_p90'].map(transform),color,true)];
 const pairEndpoints=(title,subtitle,groups,notes,source)=>{
  const s=slide(title,subtitle,notes,source), series=[];
  const labels=Object.keys(groups);
  labels.forEach((label,i)=>{
   const v=groups[label];series.push({name:label,xValues:v.map((_,j)=>i+((j%7)-3)*.025),values:v,fill:(labels.length===3?['#657580','#236B8E','#A74F3D']:['#236B8E','#A74F3D'])[i],marker:{symbol:'circle',size:6},line:{fill:'none',width:0}});
  });
  series.push({name:'中位數',xValues:labels.map((_,i)=>i),values:labels.map(k=>{const a=[...groups[k]].sort((a,b)=>a-b);return(a[14]+a[15])/2;}),fill:'#173042',marker:{symbol:'circle',size:13},line:{fill:'none',width:0}});
  const ch=scientificChart(s,'scatter',[],series,'第 18 代細胞數',{xTitle:'條件編號：對照上方圖例'});
  ch.xAxis={...ch.xAxis,min:-.2,max:labels.length-.8,majorUnit:labels.length===3?1:.2,numberFormatCode:labels.length===3?'0':'0.0'};
  side(s,'n = 30','每個點為一次試驗\n黑色大點為中位數\n\n中位差異與\n試驗間重疊一起看');
  foot(s,'固定第 18 代比較。橫向抖動僅減少遮蔽，不表示其他變數。');
 };
 {
  const s=slide('治療能降低負荷，能否持續控制？','先看總量，再追問族群組成與變異的來源',
   '這個問題承接前半部程式規則。治療降低敏感細胞存活，但抗藥細胞仍可增長。因此需用負荷、抗藥比例、多樣性與事件分布共同判讀。以下結果皆來自30次預設試驗，曲線按絕對generation對齊，每trial治療起始可能不同。', 'outputs/trajectory.csv, therapy_replicates.csv');
  line(s,'負荷下降後為何再增長？',240,75,1130,35);
  line(s,'抗藥比例升高時，哪些突變也在變化？',345,75,1130,32);
  line(s,'改變 TSG 門檻或共享環境，結果是否跟著改變？',455,75,1130,30);
  foot(s,'證據順序：治療反應 → 組成與多樣性 → 控制實驗 → 後續驗證。');
 }
 {
  const s=slide('治療後總量下降，抗藥族群再擴張','TSG two-hit、IL11 OFF，前 30 代每代 n=30',
   '先在原始細胞數上計算percentiles再轉換log10(N+1)。每組治療開始世代中位15，各自達到門檻才治療，所以不畫一條共同啟動線。總、敏感、抗藥各自中位數不必相加。813.5是每個trial最低負荷的中位數，不是中位曲線最低點。', 'outputs/trajectory.csv, therapy_replicates.csv');
  scientificChart(s,'line',r.generation,[curve('總量 P10',r.total_p10.map(x=>Math.log10(x+1)),'#A5B4BD',true),curve('總量 中位',r.total_median.map(x=>Math.log10(x+1)),'#173042'),curve('總量 P90',r.total_p90.map(x=>Math.log10(x+1)),'#A5B4BD',true),curve('敏感 中位',r.sensitive_median.map(x=>Math.log10(x+1)),'#236B8E'),curve('抗藥 中位',r.resistant_median.map(x=>Math.log10(x+1)),'#A74F3D')],'log₁₀(細胞數+1)',{max:6});
  side(s,'813.5','每次試驗最低負荷\n的中位數\n\n治療使敏感群下降\n抗藥群仍可擴張');
  foot(s,'實線為中位數，總量虛線為 P10/P90，表示試驗間變異。');
 }
 {
  const s=slide('抗藥比例迅速升高，解釋再擴張的組成','治療改變存活差異，突變率 μ 保持不變',
   'resistant_percent在存活且觀察中的trial計算。第15代中位1.1585%、第18代90.5163%、第20代99.9293%。抗藥比例增加不等同突變率增加，也不證明所有細胞同一家系。按代數保留n，這個窗口全部30個trial存活。', 'outputs/trajectory.csv');
  scientificChart(s,'line',r.generation,percentiles('',r,'#A74F3D','resistant_percent'),'抗藥細胞 (%)',{max:100});
  side(s,'99.93%','第 20 代中位比例\n\n選擇壓力使\n抗藥細胞佔優');
  foot(s,'比例依存活腫瘤定義；30 次試驗都觀察到第 31 代，之後因計算上限截尾。');
 }
 {
  const s=slide('有效基因型數顯示治療後多樣性改變','只看總負荷，會漏掉族群內部的結構變化',
   'effective_genotypes=exp(Shannon)，有相同Shannon的均勻組成需這麼多基因型。它描述基因型多樣性，不是祖先家系數。直接使用trajectory已記錄trial百分位，前30代全部trial仍觀察且存活。圖中突變持續、不同基因型生長與選擇共同影響多樣性，不能僅由負荷推斷。', 'outputs/trajectory.csv: effective_genotypes');
  scientificChart(s,'line',r.generation,percentiles('',r,'#6C5B9C','effective_genotypes'),'有效基因型數 exp(H)');
  side(s,'exp(H)','H 同時考慮\n種類與均勻程度\n\n基因型數目\n與負荷分開解讀');
  foot(s,'中位與 P10/P90，n=30。多樣性描述基因型組成，不追蹤家系。');
 }
 {
  const g=results.gene_frequency, trans=x=>Math.log10(x+.01);
  const s=slide('基因頻率：被篩選的基因與中性對照','30 次試驗以治療開始對齊，每條線是各基因突變細胞比例的中位數',
   '取代原本兩張熱圖。資料來自plot_genotype_evolution.py的gene_frequency_summary.csv：重現30個therapy seeds，以各自的治療開始為第0代對齊，只取所有試驗都有觀察的相對世代(−14到+17)。passenger mean是14個中性位點的平均，作為「沒有被篩選」的基準。治療前oncogene約4%而passenger約1.5%，突變率相同，差距來自選擇；RES在治療前與passenger同一水準，證明抗藥在治療前是中性的，治療後約4代升到100%。TSG任一擊約為passenger兩倍，因為有兩個等位基因。縱軸是log10(百分比+0.01)，0.01讓零值可以顯示。',
   'outputs/genotype_evolution/gene_frequency_summary.csv');
  scientificChart(s,'line',g.relative_generation,[curve('ONC1',g.ONC1.map(trans),'#236B8E'),curve('ONC2',g.ONC2.map(trans),'#4F8FB0'),
   curve('ONC3',g.ONC3.map(trans),'#86B5CC'),curve('RES',g.RES.map(trans),'#A74F3D'),curve('TSG1 任一擊',g['TSG1 any hit'].map(trans),'#6C5B9C'),
   curve('passenger 平均',g['passenger mean'].map(trans),'#657580',true)],'log₁₀(突變細胞 % + 0.01)',{xTitle:'相對治療開始的世代',min:-2,max:2.1,yFormat:'0.0',majorUnit:1});
  side(s,'4% vs 1.5%','治療開始時\noncogene vs passenger\n\n同突變率，差距\n來自篩選');
  foot(s,'中位數，n=30。RES 在治療前與 passenger 同一水準：抗藥性在治療前是中性的。');
 }
 {
  const hist=r.relapse_histogram;
  const s=slide('30 次試驗均復發，但事件時間仍有差異','事件定義：先低於起始治療負荷，再回到該負荷',
   '先觀察組成機制再看事件分布。relapse_generation範圍23到31，中位28，P10/P90為26/29.1。這批30/30不代表任何參數下都必然復發，更不能推論病人必然復發。未觀察事件需有截尾資料與其他分類，本批均完整觀察復發後才celllimit。', 'outputs/therapy_replicates.csv');
  scientificChart(s,'bar',hist.generation,[{name:'復發試驗數',values:hist.count,fill:'#A74F3D',valuesFormatCode:'0'}],'試驗數',{xTitle:'首次復發世代',max:10,labels:true});
  side(s,'30/30','觀察到復發\n\n事件時間中位 28\nP10/P90=26/29.1');
  foot(s,'樣本數為 30，結果限於本批教學參數與觀察窗口。');
 }
 {
  const a=results.therapy_analysis, gc=a.goldie_coldman;
  const s=slide('越早治療越可能治癒：Goldie–Coldman','治療起點 100–5,000 細胞，每個大小 200 次；理論沒有任何擬合參數',
   `analyze_therapy.py 08。橫軸log10(治療開始負荷)，縱軸治癒(滅絕)百分比。理論：抗藥突變clone以Poisson過程出現，治療前s代出生的clone最終存活機率1−f^(s)(q_on)，q_on=0.479是治療下一個抗藥家系的滅絕機率，治癒≈exp(−Λ)。只用基準細胞，所以大負荷時略高估治癒。理論50%點約${a.half_cure_size}細胞。主實驗20,000細胞遠超過這個範圍，因此30/30復發。模擬誤差以Wilson 95%區間記錄於CSV。`,
   'outputs/therapy_analysis/goldie_coldman_summary.csv');
  const x=gc.map(r=>Math.log10(r.threshold));
  const ch=scientificChart(s,'scatter',[],[
   {name:'模擬（200 次）',xValues:x,values:gc.map(r=>100*r.p_cure),fill:'#A74F3D',marker:{symbol:'circle',size:9},line:{fill:'none',width:0}},
   {name:'分枝過程理論',xValues:x,values:gc.map(r=>100*r.theory_p_cure),fill:'#173042',marker:{symbol:'diamond',size:8},line:{fill:'#173042',width:1.5}}],
   '治癒機率 (%)',{xTitle:'log₁₀(開始治療時的細胞數)',max:100});
  ch.xAxis={...ch.xAxis,min:1.9,max:3.8,majorUnit:.5,numberFormatCode:'0.0'};
  side(s,`≈ ${a.half_cure_size}`,'理論 50% 治癒\n所需的細胞數\n\n100 細胞：99.5%\n1,000 細胞：10%');
  foot(s,'治癒 = 腫瘤滅絕。治療越晚，治療前已存在的抗藥細胞越多。');
 }
 {
  const a=results.therapy_analysis, r=a.relapse_rows;
  const s=slide('復發時間在治療開始時大致已經決定','用治療開始那一代的抗藥基因型預測復發延遲，30 次試驗',
   `analyze_therapy.py 09。基準預測ln(N/R0)/ln 1.371，假設所有抗藥細胞相同，平均誤差${a.relapse.mae_baseline.toFixed(1)}代且偏晚。組成預測讓每種抗藥基因型以自己的治療中平均子代數成長(帶oncogene ×1.62)，敏感細胞×0.143衰減；向上取整後30次中${a.relapse.exact_genotype}次完全命中，平均誤差${a.relapse.mae_genotype.toFixed(1)}代。R0越多越早復發(Spearman ${a.relapse.spearman_obs_vs_R0.toFixed(2)})，oncogene比例越高越早(${a.relapse.spearman_obs_vs_onc.toFixed(2)})。治療中成長率在ρ*=${(100*a.relapse.rho_star).toFixed(0)}%穿過零，即谷底。`,
   'outputs/therapy_analysis/relapse_time_prediction.csv');
  const ch=scientificChart(s,'scatter',[],[
   {name:'基準預測',xValues:r.map(v=>v.predicted_delay_baseline),values:r.map(v=>v.observed_relapse_delay),fill:'#A5B4BD',marker:{symbol:'circle',size:7},line:{fill:'none',width:0}},
   {name:'考慮基因型組成',xValues:r.map(v=>v.predicted_delay_genotype),values:r.map(v=>v.observed_relapse_delay),fill:'#236B8E',marker:{symbol:'circle',size:9},line:{fill:'none',width:0}}],
   '觀察到的復發延遲（代）',{xTitle:'預測的復發延遲（治療開始後的代數）',min:7,max:18});
  ch.xAxis={...ch.xAxis,min:7,max:18,majorUnit:2,numberFormatCode:'0'};
  side(s,`${a.relapse.exact_genotype}/30`,`完全命中的試驗\n\n平均誤差\n${a.relapse.mae_genotype.toFixed(1)} 代（基準 ${a.relapse.mae_baseline.toFixed(1)}）\n\n谷底：ρ* ≈ ${(100*a.relapse.rho_star).toFixed(0)}%`);
  foot(s,'帶 oncogene 的抗藥細胞長得快，是基準預測偏晚的原因。');
 }
 {
  const a=results.therapy_analysis.sensitivity, keep=[0,2,4,5,6,7,8];
  const s=slide('殺傷率要壓過最強的抗藥細胞才能治癒','抗藥細胞殺傷率 k_R × 突變率 μ 的治癒機率，每格 30 次',
   `analyze_therapy.py 10。表中數字是治癒(滅絕)百分比；預設k_R=0.04、μ=0.001。單純抗藥細胞在k_R=${a.boundaries['baseline resistant'].toFixed(2)}停止成長，但帶oncogene的抗藥細胞要到${a.boundaries['+ oncogene'].toFixed(2)}，再加TSG失活是${a.boundaries['+ oncogene + TSG loss'].toFixed(2)}。突變率越高，治療前越可能已有這種雙突變細胞，治癒邊界越往右移。完整熱圖與復發延遲見outputs/therapy_analysis/10_parameter_sensitivity.png。`,
   'outputs/therapy_analysis/parameter_sensitivity_trials.csv');
  const rows=[['μ ＼ k_R',...keep.map(i=>String(a.kills[i]))]];
  for(let j=a.mutation_rates.length-1;j>=0;j--)rows.push([String(a.mutation_rates[j]),...keep.map(i=>Math.round(100*a.p_cure[j][i])+'%')]);
  table(s,rows,70,215,880,330,[160,...keep.map(()=>102)],23);
  side(s,'0.30 → 0.41','單純抗藥 → 帶\noncogene 的抗藥\n細胞停止成長\n的殺傷率');
  foot(s,'治癒邊界由最強的抗藥基因型決定，不是一般的抗藥細胞。');
 }
 {
  const a=results.therapy_analysis.adaptive;
  const s=slide('適應性治療只有在細胞互相競爭時才有用','間歇用藥 vs 持續用藥的惡化時間中位數，每種設定 30 對試驗',
   `analyze_therapy.py 11。間歇治療：負荷降到治療起點一半就暫停，回到20,000再恢復。惡化=負荷曾低於20,000後，在用藥中仍≥24,000。兩種策略用同一組seed。沒有競爭時敏感細胞不會妨礙抗藥細胞，暫停只是讓抗藥細胞少挨殺傷，所以間歇較差。加入承載量競爭(K=30,000，每個細胞額外死亡0.21·N/K)後間歇延後惡化，配對比較中幾乎每次都較長，且只用約一半藥量。代價(抗藥細胞死亡+0.05)讓差距更明顯。競爭形式與參數都是教學假設。`,
   'outputs/therapy_analysis/adaptive_therapy_summary.csv');
  const labels=['無競爭','無競爭＋代價','競爭','競爭＋代價'];
  scientificChart(s,'bar',labels,[{name:'持續用藥',values:a.map(r=>r.continuous_median_ttp),fill:'#236B8E',valuesFormatCode:'0.0'},
   {name:'間歇（適應性）',values:a.map(r=>r.adaptive_median_ttp),fill:'#A74F3D',valuesFormatCode:'0.0'}],'惡化時間中位數（代）',{xTitle:'模型設定',labels:true,max:24});
  side(s,'×1.14','有競爭時\n間歇用藥延後惡化\n\n無競爭時\n反而較差（×0.85）');
  foot(s,'主模型沒有競爭，因此在主模型中持續用藥較好。');
 }
 {
  const c=results.no_tsg_advantage,a=results.hits['TSG one-hit'],b=results.hits['TSG two-hit'];
  const s=slide('控制實驗：TSG 功能優勢改變生長','保留 oncogene 與突變設定，比較無 TSG 優勢、one-hit 和 two-hit',
   '對照組保持相同基因組及oncogene規則，只把TSG deathreduction和triplebonus設0。onehit/twohit只改同一TSG門檻。新重現30個對照seed並核對comparisons.csv終點。三條中位線按18代固定時間，不治療，避免治療選擇乾擾。此圖對應02_scenario_comparisons左面板。', 'outputs/comparisons.csv, one_hit_vs_two_hit.csv');
  scientificChart(s,'line',c.generation,[curve('無 TSG 優勢',c.total_median,'#657580'),curve('one-hit',a.total_median,'#236B8E'),curve('two-hit',b.total_median,'#A74F3D')],'細胞數');
  side(s,'控制變因','相同突變機會\n相同 oncogene 規則\n\n差異來自 TSG\n功能與達標門檻');
  foot(s,'每組 30 次，無治療，固定 18 代。下一頁保留試驗變異範圍。');
 }
 {
  const a=results.hits['TSG one-hit'],b=results.hits['TSG two-hit'];
  const s=slide('one-hit 組的第 18 代中位負荷較高','比較同一時間點，並呈現兩組試驗變異的重疊',
   'onehit中位94617.5，twohit73879，兩組中位比1.2807，即高28.1%。不是每個onehittrial都大於每個twohittrial。這是功能門檻改變下的教學模型結果，無pvalue聲稱。', 'outputs/one_hit_vs_two_hit.csv');
  scientificChart(s,'line',a.generation,[...percentiles('one-hit',a,'#236B8E'),...percentiles('two-hit',b,'#A74F3D')],'細胞數',{max:160000});
  side(s,'+28.1%','第 18 代中位負荷比\n\n94,617.5 / 73,879\n\n變異範圍仍重疊');
  foot(s,'各組 n=30。實線中位，虛線 P10/P90，非信賴區間。');
 }
 {
  const a=results.hits['TSG one-hit'],b=results.hits['TSG two-hit'],c0=results.hits['No TSG advantage'];
  const trans=x=>Math.log10(100*x+.01);
  const s=slide('TSG 達標比例揭示門檻差異的機制','這裡計算功能喪失比例，與基因熱圖的 any hit 分開',
   'onehit只要任一等位基因即可達標，twohit同基因兩等位基因，所以達標比例可能跨量級。為保留0並可顯示稀少twohit，顯式轉換log10(p_percent+.01)；0%=-2，.99%=0，9.99%=1。原始fig03面板B用symlog；這裡改用明確注明的平滑對數表示而非假裝沿用同一軸。資料用suppressor_fraction。', 'outputs/one_hit_vs_two_hit.csv: suppressor_fraction');
  scientificChart(s,'line',a.generation,[curve('one-hit 中位',a.suppressor_fraction_median.map(trans),'#236B8E'),curve('two-hit 中位',b.suppressor_fraction_median.map(trans),'#A74F3D'),curve('無 TSG 優勢 中位',c0.suppressor_fraction_median.map(trans),'#657580',true)],'log₁₀(TSG 達標 %+0.01)',{min:-2,max:2,yFormat:'0.0',majorUnit:1});
  side(s,'同基因兩擊','two-hit 達標需\n累積第二等位基因\n\n0% 對應 −2\n0.99% 對應 0');
  foot(s,'各組 n=30，中位曲線。加 0.01 後轉換，避免零值無法取對數。');
 }
 pairEndpoints('TSG 終點分布：two-hit 與對照組幾乎相同','條件 0：無 TSG 優勢　1：one-hit　2：two-hit，每組 n=30',results.endpoints.hits,
  '這是原始fig03面板C的固定18代endpoint散點。所有trial逐點保留，加上中位數。兩組不同seeds，並非配對點。可看到分布重疊，避免把中位28.1%差異當成每一個trial都相同。', 'outputs/one_hit_vs_two_hit.csv: generation 18');

 {
  const c=results.il11_curve;
  const s=slide('IL11 機制假設：少量生產者即可形成共享效益','先看模型規則，再看 ON/OFF 模擬是否符合預期',
   '這張是原始fig05低比例放大面板。曲線來自指定公式而不是隨機模擬或實驗擬合。f=.02時b=.03，f=.1時b=.05。bmax=.06，不會超過6個百分點。半飽和參數是教學假設，可能影響結果。所有細胞均受益，未考慮分泌成本或空間限制。', 'outputs/05_il11_benefit_curve.svg, Config');
  scientificChart(s,'line',c.fraction_percent.filter((_,i)=>i%5===0),[curve('共享死亡降低',c.benefit_percentage_points.filter((_,i)=>i%5===0),'#236B8E')],'死亡率降低 (百分點)',{xTitle:'IL11⁺ 生產者比例 (%)',max:6});
  side(s,'f = 2%','b = 3 個百分點\n達最大效益的一半\n\n實際生產者 0–2%\n第 18 代中位 1.6%\n正落在最陡區段');
  foot(s,'顯示 0–10% 生產者區間。上限為 6 個百分點，完整範圍漸近此上限。');
 }
 {
  const a=results.il11['IL11 OFF'],b=results.il11['IL11 ON'];
  const s=slide('IL11 ON 組的中位負荷提高 52.5%','相同位點與突變機會，啟用共享死亡降低效益',
   'OFF中位72455.5、ON110469.5，比值1.52465。每rep同seed，但族群分歧後隨機數消耗不同，不是逐細胞配對。控制相同20loci避免突變機會變化。未治療18代。52.5%是兩組中位比，不是pairedtrial提升平均。', 'outputs/il11_comparison.csv');
  scientificChart(s,'line',a.generation,[...percentiles('OFF',a,'#236B8E'),...percentiles('ON',b,'#A74F3D')],'細胞數',{max:180000});
  side(s,'+52.5%','第 18 代中位負荷比\n\n共享存活效益\n使整個族群受益');
  foot(s,'各組 n=30，實線中位、虛線 P10/P90。其他參數固定，不施加治療。');
 }
 {
  const a=results.il11['IL11 OFF'],b=results.il11['IL11 ON'];
  const s=slide('生產者仍是少數，共享效益已可影響全族群','把生產者比例與負荷差異連起來解釋',
   '對應fig04面板B。il11_percent是生產者比例，中位ON第18代1.603%，OFF位點中性但仍可能突變。ON所算shareddeathreduction中位2.6694百分點。比例中位與benefit中位因單調轉換一致但為每trial計算，不能把生產者數少誤解為只有少數細胞受益。', 'outputs/il11_comparison.csv: il11_percent, shared_death_reduction');
  scientificChart(s,'line',a.generation,[...percentiles('OFF',a,'#236B8E','il11_percent'),...percentiles('ON',b,'#A74F3D','il11_percent')],'IL11 位點突變細胞 (%)',{yFormat:'0.0'});
  side(s,'1.60%','ON 第 18 代中位\n\n對應死亡率降低\n約 2.67 個百分點\n\n非生產者同樣受益');
  foot(s,'OFF 為同一中性位點的突變比例。兩組保留相同位點與機會。');
 }
 pairEndpoints('IL11 終點散點：共享效益與隨機變異並存','條件 0：OFF　條件 1：ON，固定第 18 代',results.endpoints.il11,
  '原始fig04面板C。同seed保持experiment可重現，但不能當成逐細胞配對。散點呈現原始30個endpoint，黑點為中位。看差異也看重疊，不聲稱實驗資料校準或統計顯著性。', 'outputs/il11_comparison.csv: generation 18');
 {
  const s=slide('高死亡率組的等大小樣本較多樣','首次達標後無放回取 5,000 細胞，每組 30 次',
   '原始fig02右面板。d=.2,.35,.45均30/30到達取樣目標。Shannon中位1.06236,1.80123,2.83497；達標generation中位8,14,26.5。較慢達標可累積更多mutation，因此不能單獨認定death提高多樣性。取樣控制的是樣本數，沒有控制演化時間。', 'outputs/comparisons.csv: death_rate');
  const series=results.death_rate.map((d,i)=>({name:'d='+d.death,xValues:d.shannon.map((_,j)=>d.death+((j%7)-3)*.002),values:d.shannon,fill:['#236B8E','#A74F3D','#6C5B9C'][i],marker:{symbol:'circle',size:6},line:{fill:'none',width:0}}));
  series.push({name:'中位數',xValues:results.death_rate.map(d=>d.death),values:results.death_rate.map(d=>d.median_shannon),fill:'#173042',marker:{symbol:'circle',size:13},line:{fill:'none',width:0}});
  const ch=scientificChart(s,'scatter',[],series,'Shannon H',{xTitle:'自然死亡機率 d',max:4});ch.xAxis={...ch.xAxis,min:.15,max:.5,majorUnit:.05,numberFormatCode:'0.00'};
  side(s,'1.06 / 1.80\n2.83','三組 H 中位數\n\n達標世代中位數\n8、14、26.5');
  foot(s,'每點為一次試驗，黑色大點為中位數。比較等大小樣本，達標時間仍不同。');
 }
 {
  const s=slide('結果支持的機制與下一步實驗','每一個解釋，都對應可改動的程式參數和可量測指標',
   '抗藥殺傷率掃描用replace修改resistanttherapy kill，其餘保持固定。IL11掃描bmax與K，加入producercost或spatialrestriction。多樣性固定代數及固定樣本大小分開比較並保留未達標trial。先定義指標再新增experiment，存CSV與metadata及tests。需生物資料校準後再談實際適用。');
  table(s,[['觀察','程式中的可能機制','下一步與指標'],['抗藥佔優與復發','kill_R 要壓過最強抗藥型（已掃描）','加入免疫或空間結構\n比較治癒邊界是否移動'],['IL11 ON 負荷較高','b(f) 讓非生產者也受益','掃描 b_max、K，加分泌成本\n負荷差異、生產者比例'],['高死亡組 H 較高','周轉與達標時間共同影響','固定世代和負荷分別比較\nH、達標率與代數']],70,215,1140,355,[280,350,510],25);
  foot(s,'這些是指定參數下的模擬觀察。後續需敏感度分析與資料校準。');
 }
}
