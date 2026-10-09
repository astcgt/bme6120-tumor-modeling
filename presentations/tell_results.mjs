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
   const v=groups[label];series.push({name:label,xValues:v.map((_,j)=>i+((j%7)-3)*.025),values:v,fill:['#236B8E','#A74F3D'][i],marker:{symbol:'circle',size:6},line:{fill:'none',width:0}});
  });
  series.push({name:'中位數',xValues:[0,1],values:labels.map(k=>{const a=[...groups[k]].sort((a,b)=>a-b);return(a[14]+a[15])/2;}),fill:'#173042',marker:{symbol:'circle',size:13},line:{fill:'none',width:0}});
  const ch=scientificChart(s,'scatter',[],series,'第 18 代細胞數',{xTitle:'條件：0 與 1 對照上方圖例'});
  ch.xAxis={...ch.xAxis,min:-.2,max:1.2,majorUnit:.2,numberFormatCode:'0.0'};
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
  foot(s,'比例依存活腫瘤定義。P10/P90 是變異範圍，並非信賴區間。');
 }
 {
  const s=slide('有效基因型數顯示治療後多樣性改變','只看總負荷，會漏掉族群內部的結構變化',
   'effective_genotypes=exp(Shannon)，有相同Shannon的均勻組成需這麼多基因型。它描述基因型多樣性，不是祖先家系數。直接使用trajectory已記錄trial百分位，前30代全部trial仍觀察且存活。圖中突變持續、不同基因型生長與選擇共同影響多樣性，不能僅由負荷推斷。', 'outputs/trajectory.csv: effective_genotypes');
  scientificChart(s,'line',r.generation,percentiles('',r,'#6C5B9C','effective_genotypes'),'有效基因型數 exp(H)');
  side(s,'exp(H)','H 同時考慮\n種類與均勻程度\n\n基因型數目\n與負荷分開解讀');
  foot(s,'中位與 P10/P90，n=30。多樣性描述基因型組成，不追蹤家系。');
 }
 {
  const a=results.availability;
  const s=slide('後期可用試驗數下降，是觀察截尾','完整 0–40 代資料讓停止規則的影響可見',
   '這是原始dashboard面板C。celllimit停止後不能再列為觀察，曲線中的observed_n下降。不把它誤解成治療成功。所有30次試驗最終celllimit，沒有extinction。therapy_started計算每代已達啟動門檻的trial累計。主曲線與heatmap限制0-30，是因為該窗口全部30仍存活觀察。', 'outputs/trajectory.csv, run_metadata.json');
  scientificChart(s,'line',a.generation,[curve('有效負荷 n',a.observed_n,'#236B8E'),curve('存活且觀察 n',a.alive_n,'#A74F3D'),curve('已啟動治療',a.therapy_started,'#6C5B9C',true)],'試驗數',{max:32});
  side(s,'n 隨代變','達計算上限後\n後續數值缺失\n\n解讀後期曲線\n先確認有效 n');
  foot(s,'停止觀察與滅絕有不同意義。有效負荷與存活 n 本批次重疊。');
 }
 const heatmap=(mode,title,subtitle)=>{
  const h=results.gene_heatmaps,values=h[mode],max=mode==='mean'?1:.5;
  const s=slide(title,subtitle,
   `重現30個therapy seeds並逐代核對therapy_trajectories.csv的total/resistant共1119列。基因頻率是各trial中該位點任一等位基因突變的比例。TSG顯示any hit，不等於two-hit功能喪失。${mode==='mean'?'每格為30次試驗的平均頻率。':'每格為30次試驗頻率的樣本SD，ddof=1。'}採0-30代的完整共同窗口。色階${max}。頻率之間可能因共現而關聯，這張圖不提供因果關係或家系識別。`,
   'outputs/01_tumor_dashboard panels E/F; regenerated populations from recorded seeds');
  const x0=218,y0=216,cw=23.4,rh=17.3;
  for(let i=0;i<20;i++){
   text(s,h.gene_labels[i],72,y0+i*rh,140,rh,13.5,i===3?'#A74F3D':'#173042',i<6);
   values[i].forEach((v,j)=>{
    const a=Math.max(0,Math.min(1,v/max));
    const rgb=[248,251,253].map((c,k)=>Math.round(c+([35,107,142][k]-c)*a));
    s.shapes.add({geometry:'rect',position:{left:x0+j*cw,top:y0+i*rh,width:cw+.1,height:rh+.1},fill:'#'+rgb.map(v=>v.toString(16).padStart(2,'0')).join(''),line:{fill:'none',width:0}});
   });
  }
  for(let g=0;g<=30;g+=5)text(s,String(g),x0+g*cw-5,568,48,25,18);
  text(s,'世代',545,595,100,28,20);
  side(s,mode==='mean'?'平均頻率':'試驗間 SD',mode==='mean'?'抗藥位點逐漸固定\n其他位點一同變化\n\nTSG 顏色表示\n任一等位基因突變':'低平均頻率也可\n有較大試驗差異\n\nSD 描述隨機變異\n不代表測量誤差');
  for(let j=0;j<20;j++){
   const a=j/19,rgb=[248,251,253].map((c,k)=>Math.round(c+([35,107,142][k]-c)*a));
   s.shapes.add({geometry:'rect',position:{left:980+j*10.5,top:375,width:11,height:16},fill:'#'+rgb.map(v=>v.toString(16).padStart(2,'0')).join(''),line:{fill:'none',width:0}});
  }
  text(s,'0',980,394,35,30,18);
  text(s,String(max),1158,394,32,30,18,'#173042',false,'right');
  foot(s,`每格 n=30。${mode==='mean'?'平均':'樣本 SD'}色階 0–${max}，TSG 為 any hit，並非功能喪失比例。`);
 };
 heatmap('mean','基因頻率熱圖：哪些位點隨治療改變？','20 個基因座，前 30 代每格為 30 次試驗平均');
 heatmap('sd','基因頻率 SD：組成並非每次都相同','相同參數下仍有試驗間變異，需與平均熱圖一起讀');
 {
  const hist=r.relapse_histogram;
  const s=slide('30 次試驗均復發，但事件時間仍有差異','事件定義：先低於起始治療負荷，再回到該負荷',
   '先觀察組成機制再看事件分布。relapse_generation範圍23到31，中位28，P10/P90為26/29.1。這批30/30不代表任何參數下都必然復發，更不能推論病人必然復發。未觀察事件需有截尾資料與其他分類，本批均完整觀察復發後才celllimit。', 'outputs/therapy_replicates.csv');
  scientificChart(s,'bar',hist.generation,[{name:'復發試驗數',values:hist.count,fill:'#A74F3D',valuesFormatCode:'0'}],'試驗數',{xTitle:'首次復發世代',max:10,labels:true});
  side(s,'30/30','觀察到復發\n\n事件時間中位 28\nP10/P90=26/29.1');
  foot(s,'樣本數為 30，結果限於本批教學參數與觀察窗口。');
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
  const a=results.hits['TSG one-hit'],b=results.hits['TSG two-hit'];
  const trans=x=>Math.log10(100*x+.01);
  const s=slide('TSG 達標比例揭示門檻差異的機制','這裡計算功能喪失比例，與基因熱圖的 any hit 分開',
   'onehit只要任一等位基因即可達標，twohit同基因兩等位基因，所以達標比例可能跨量級。為保留0並可顯示稀少twohit，顯式轉換log10(p_percent+.01)；0%=-2，.99%=0，9.99%=1。原始fig03面板B用symlog；這裡改用明確注明的平滑對數表示而非假裝沿用同一軸。資料用suppressor_fraction。', 'outputs/one_hit_vs_two_hit.csv: suppressor_fraction');
  scientificChart(s,'line',a.generation,[curve('one-hit 中位',a.suppressor_fraction_median.map(trans),'#236B8E'),curve('two-hit 中位',b.suppressor_fraction_median.map(trans),'#A74F3D')],'log₁₀(TSG 達標 %+0.01)',{min:-2,max:2,yFormat:'0.0',majorUnit:1});
  side(s,'同基因兩擊','two-hit 達標需\n累積第二等位基因\n\n0% 對應 −2\n0.99% 對應 0');
  foot(s,'各組 n=30，中位曲線。加 0.01 後轉換，避免零值無法取對數。');
 }
 pairEndpoints('TSG 終點分布：兩組試驗仍有重疊','條件 0：one-hit　條件 1：two-hit，每組 n=30',results.endpoints.hits,
  '這是原始fig03面板C的固定18代endpoint散點。所有trial逐點保留，加上中位數。兩組不同seeds，並非配對點。可看到分布重疊，避免把中位28.1%差異當成每一個trial都相同。', 'outputs/one_hit_vs_two_hit.csv: generation 18');

 {
  const x=Array.from({length:21},(_,i)=>i*5),y=x.map(f=>6*(f/100)/(.02+f/100));
  const s=slide('IL11 飽和曲線：共享效益有明確上限','完整 0–100% 範圍，下一頁放大實際較常見的低比例區域',
   '對應原始fig05全範圍面板。x是生產者百分比，y是死亡降低百分點，公式bmaxf/(K+f)指定上限6個百分點。f=100%時5.882百分點，接近但不等於6。這是模型假設的視覺化，不能當成模擬或生物實驗證據。低比例敏感度由下一頁放大。', 'outputs/05_il11_benefit_curve.svg: full-range panel');
  scientificChart(s,'line',x,[curve('共享死亡降低',y,'#236B8E')],'死亡率降低 (百分點)',{xTitle:'IL11⁺ 生產者比例 (%)',max:6});
  side(s,'上限 6','高比例時效益趨緩\n\nb_max 決定上限\nK 決定低比例敏感度');
  foot(s,'本圖來自指定飽和函數。生產者比例 100% 時降低量約 5.88 個百分點。');
 }
 {
  const c=results.il11_curve;
  const s=slide('IL11 機制假設：少量生產者即可形成共享效益','先看模型規則，再看 ON/OFF 模擬是否符合預期',
   '這張是原始fig05低比例放大面板。曲線來自指定公式而不是隨機模擬或實驗擬合。f=.02時b=.03，f=.1時b=.05。bmax=.06，不會超過6個百分點。半飽和參數是教學假設，可能影響結果。所有細胞均受益，未考慮分泌成本或空間限制。', 'outputs/05_il11_benefit_curve.svg, Config');
  scientificChart(s,'line',c.fraction_percent.filter((_,i)=>i%5===0),[curve('共享死亡降低',c.benefit_percentage_points.filter((_,i)=>i%5===0),'#236B8E')],'死亡率降低 (百分點)',{xTitle:'IL11⁺ 生產者比例 (%)',max:6});
  side(s,'f = 2%','b = 3 個百分點\n達最大效益的一半\n\n這是指定公式\n不是觀察資料的擬合');
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
  table(s,[['觀察','程式中的可能機制','下一步與指標'],['抗藥佔優與復發','kill_S 與 kill_R 的差異','掃描 resistant_therapy_kill\n復發率、時間、滅絕率'],['IL11 ON 負荷較高','b(f) 讓非生產者也受益','掃描 b_max、K，加分泌成本\n負荷差異、生產者比例'],['高死亡組 H 較高','周轉與達標時間共同影響','固定世代和負荷分別比較\nH、達標率與代數']],70,215,1140,355,[280,350,510],25);
  foot(s,'這些是指定參數下的模擬觀察。後續需敏感度分析與資料校準。');
 }
}
