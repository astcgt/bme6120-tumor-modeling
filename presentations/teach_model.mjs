/** Code-first lessons. Snippets are selected verbatim from modeling_cancer.py. */
export function teachModel({slide,text,line,foot,node,connect,table,math,code,snip,P,pi2,scientificChart,curve,results}){
 const lesson=(title,subtitle,range,explanations,formula,notes)=>{
  const s=slide(title,subtitle,notes,`modeling_cancer.py:${range.map(r=>r.join('-')).join(', ')}`);
  code(s,snip(range),70,215,650,300);
  explanations.forEach((p,i)=>text(s,p,770,219+i*97,435,92,24));
  if(formula)math(s,formula,70,formula.includes('\\frac')||formula.includes('\\sum')?520:550,1135,formula.includes('\\frac')||formula.includes('\\sum')?95:65);
  foot(s,'原始碼節錄，長行採視覺折行。省略註解與非核心行，詳細說明見備註。');
  return s;
 };
 {
  const s=slide('腫瘤模擬程式的設計與數學','BME 6120　程式實作、狀態轉移與模擬結果',
   '本簡報以實作為主線。先建立資料與參數的語意，再追蹤simulate的一代更新，接著理解突變抽樣與Markov核，最後追溯結果、統計處理與測試。所有數學在對應概念頁直接說明，沒有附錄。結果沿用已核對的30次模擬資料。');
  line(s,'一代更新：存活、分裂，再讓子代突變',245,75,1130,36);
  math(s,String.raw`S_t=(\mathbf{N}_t,A_t),\qquad S_{t+1}\mid S_t=s\sim P(s,\cdot)`,70,365,1135,90);
  line(s,'閱讀順序：資料表示 → 更新規則 → 統計與證據',510,75,1130,29);
  foot(s,'程式以基因型計數模擬同步世代，參數為教學假設。');
 }

 {
  const r=results.therapy;
  const s=slide('負荷下降之後，為何又回升？','主治療情境的 30 次試驗，先提出需要解釋的現象',
   '開場圖只看總量，先提出下降後回升的研究問題。資料是TSG twohit、IL11off，30次trial，前30代完整觀察。曲線是log10總量加1之後的中位與百分位；治療各trial各自達門檻後開始。本圖本身不能回答是什麼基因型再擴張，需由接下來的程式規則與組成圖解釋。', 'outputs/trajectory.csv, run_metadata.json');
  scientificChart(s,'line',r.generation,[curve('P10',r.total_p10.map(x=>Math.log10(x+1)),'#A5B4BD',true),curve('總量中位',r.total_median.map(x=>Math.log10(x+1)),'#236B8E'),curve('P90',r.total_p90.map(x=>Math.log10(x+1)),'#A5B4BD',true)],'log₁₀(細胞數+1)',{max:6});
  text(s,'負荷下降\n組成如何改變？',970,270,240,135,30,'#236B8E',true);
  text(s,'用程式規則\n建立機制解釋\n再回到結果驗證',970,455,240,135,25);
  foot(s,'先看現象，再讀實作與數學，最後用完整圖組檢查解釋。');
 }
 {
  const s=slide('程式架構：核心模擬與實驗分析','每個函數解決一個明確問題',
   'Config集中參數，Result承載單次試驗。driver_count、suppressor_loss把基因型轉為表型，growth_probabilities將表型轉為機率，il11_environment計算共享環境。simulate執行同步更新，mutate_offspring重新分類子代。summarize只讀族群統計。main安排多次試驗，plot_therapy_ensemble、run_comparisons、run_il11_comparison產生圖與CSV，therapy_outcome從history識別復發。',
   'modeling_cancer.py:32-255, 258-279, 639-762');
  const a=node(s,'Config\n參數與驗證',70,218,235,90),b=node(s,'simulate()\n一條隨機路徑',390,218,290,90),c=node(s,'Result\n摘要與完整快照',790,218,360,90);connect(s,a,b);connect(s,b,c);
  text(s,'表型與機率\ndriver_count / suppressor_loss\ngrowth_probabilities / il11_environment',70,360,635,130,25);
  text(s,'多次試驗與結果\nmain / plot_therapy_ensemble\ntherapy_outcome / run_comparisons',755,360,455,130,25);
  foot(s,'核心模擬決定資料如何生成，分析函數決定如何解讀資料。');
 }
 lesson('Config：把假設變成可檢查的參數','設定集中管理，避免各函數使用不同規則',[[32,36],[49,54],[60,60]],
  ['frozen=True 保護本次設定\nreplace() 建立新實驗條件','μ 是每個未突變位元的機率\nq 是存活後產生三子代的機率','max_cells 是計算停止門檻\n不代表腫瘤承載量'],
  String.raw`\mu=0.001,\quad d_0=0.30,\quad q_0=0.04,\quad N_0=100`,
  'Config使用不可變dataclass。n_genes是基因座數，而非總位元數。driver三個位點、resistance一個位點、TSG兩個位點和IL11一個保留位點互不重疊。validate檢查每個機率，以及三子代機率加總不超過1。更換設定使用dataclasses.replace，不修改原config。這些是模型假設，不是估計的病人參數。');
 lesson('population：稀疏字典表示整個腫瘤','字典鍵是基因型整數，值是具有該基因型的細胞數',[[233,235],[242,246]],
  ['{0: 100} 表示 100 個全未突變細胞','history 保存逐代摘要\npopulations 保存字典副本','只儲存實際出現的基因型\n不必配置 2²² 個空欄位'],
  String.raw`\mathbf{N}_t=(N_{g,t})_g,\qquad T_t=\sum_gN_{g,t}`,
  '例如population={0:90,8:10}表示90個未突變細胞與10個抗藥細胞。整數8的二進位第3位為1，因此具有resistance_gene=3。初始基因型0並不是第0個細胞，而是所有位元為0的基因型。字典副本避免後續更新污染已保存的快照。這是基因型計數，不追蹤每個細胞的家系。');
 lesson('Bitmask：用位元檢查基因是否突變','1 << g 選取第 g 位，AND 檢查，OR 設定',[[120,123]],
  ['第 3 位的遮罩為 1 << 3 = 8\n抗藥細胞滿足 genotype & 8','例：9 = 00001001₂\n第 0 與第 3 位已突變','driver_count 計數 oncogene\n抗藥與 driver 表型分開判定'],
  String.raw`I_R(g)=\mathbf{1}\{\mathrm{bit}_3(g)=1\},\qquad 9=2^0+2^3`,
  'driver_count針對config.driver_genes中的位點把布林結果加總。bool(genotype & mask)判定該位元是否非零。1<<index生成單位元遮罩。位元OR會保留已有突變，因此不可逆突變可以自然表示。抗藥突變本身不增加growth_probabilities中的生長優勢。');
 lesson('TSG：two-hit 必須擊中同一基因的兩個等位基因','20 個基因座，兩個 TSG 各增加第二等位基因位元',[[125,129],[132,136]],
  ['TSG gene 4：位元 4 與 20\nTSG gene 5：位元 5 與 21','two-hit：同一基因有兩擊\n不同基因各一擊尚未達標','any() 表示任一 TSG 達標\n每類優勢只加入一次'],
  String.raw`I_T(g)=\mathbf{1}\{\exists j:\ a_j(g)\geq h\},\qquad h\in\{1,2\}`,
  'second=config.n_genes+config.suppressor_genes.index(gene)，把TSG的第二等位基因放在20和21。suppressor_alleles把兩個位元各轉成0/1再相加。例(1<<4)|(1<<5)對two-hit為False，而(1<<4)|(1<<20)為True。one-hit條件h=1表示任一等位基因失活即可，不改變基因組大小。總位元L=20+2=22。');
 lesson('growth_probabilities：表型如何改變生長','Oncogene 與 TSG 優勢可相加，IL11 降低共享死亡率',[[140,144],[148,150]],
  ['I_O：oncogene 是否達門檻\nI_T：任一 TSG 是否達門檻','death 自然死亡率至少為 0\ntriple 是存活條件下的機率','基準 (d,q)=(0.30,0.04)\n兩類達標時 (0.14,0.28)'],
  String.raw`d_g=\max\{0,0.30-0.08I_O(g)-0.08I_T(g)-b\},\quad q_g=0.04+0.12I_O(g)+0.12I_T(g)`,
  'oncogene和suppressor是布林值，Python乘法將它們當0或1。driver達門檻降低死亡0.08並增加三子代機率0.12；TSG達門檻另加相同效果。不是每多一個driver就反覆加成。driver_hits預設1，因此任一oncogene可達標。基準均值1.428，單類優勢均值0.78*2.16=1.6848，兩類優勢0.86*2.28=1.9608。');
 lesson('il11_environment：當前族群形成下一代共享環境','本代生產者比例決定 b，所有基因型共同受益',[[156,162]],
  ['f = IL11⁺ 細胞數 / 總細胞數\n關閉或滅絕時回傳 (0,0)','f=2% 時 b=0.03\n基準死亡率由 30% 降至 27%','先計算環境，再更新子代\n本代新突變從下次轉移起作用'],
  String.raw`f_t=\frac{\sum_g N_{g,t}I_{IL11}(g)}{T_t},\qquad b_t=\frac{0.06f_t}{0.02+f_t}`,
  '均勻混合假設讓非生產者也獲得相同死亡降低。半飽和參數K=0.02，當f=K時b=bmax/2。il11_environment在每一代繁殖前呼叫一次。新子代突變在最後一步發生，所以新生產者不會回頭改變這一代已執行的存活抽樣。開關OFF保留第6位作passenger，以保持相同突變機會。');
 {
  const s=slide('simulate：一代更新的完整順序','先記錄與檢查，再繁殖與突變，最後更換族群',
   'simulate使用同步世代。記錄t代，檢查停止條件，依當前總量更新治療開關，計算IL11，再對所有父基因型存活與分裂，最後讓全部子代突變。新一代完全取代父代。把這個順序讀清楚可以避免把新突變誤當成本次存活的原因。', 'modeling_cancer.py:214-255');
  const labels=['記錄 Sₜ\n摘要與快照','停止檢查\n決定是否推進','治療與 IL11\n本次環境','父細胞抽樣\n存活與分裂','子代突變\n合併與替換'];
  const ns=labels.map((a,i)=>node(s,a,65+i*245,235,215,125));ns.slice(1).forEach((a,i)=>connect(s,ns[i],a));
  math(s,String.raw`S_t\longmapsto (\alpha_t,b_t)\longmapsto (U_g,V_g,B_g)_g\longmapsto\mathbf{N}_{t+1}`,70,420,1140,90);
  foot(s,'第 0 代也會記錄。世代更新使用全新 offspring，不把父代再加回去。');
 }
 lesson('停止條件：記錄最後一代，保留停止原因','優先序為滅絕、取樣目標、計算上限、觀察代數',[[247,257]],
  ['滅絕：族群確實已變成 0\n之後細胞數可延續為 0','達 max_cells：停止觀察\n該代完整保留，之後是缺失','target_size 用於等大小比較\n允許跨門檻時 overshoot'],
  String.raw`\tau=\inf\{t:T_t=0\ \mathrm{or}\ T_t\geq T_{\rm target}\ \mathrm{or}\ T_t\geq T_{\max}\ \mathrm{or}\ t=G\}`,
  '停止檢查在治療開關更新之前，已停止的族群不再繁殖。max_cells不是硬截斷細胞數，若上一代更新產生超過上限的細胞，程式保留完整該代再停止。計算上限不是生物吸收狀態。執行鏈若需要時間齊次狀態，可擴充(N,A,t,C)，其中C保存執行或停止模式。');
 lesson('治療開關：首次達門檻後持續作用','therapy_generation 同時保存啟動代數與控制狀態',[[258,261],[281,284]],
  ['None 表示尚未啟動治療\n達 20,000 時記錄當前代數','存活機率先乘自然存活\n再乘治療下的存活機率','敏感 k=0.90，抗藥 k=0.04\n治療改變存活，不改變 μ'],
  String.raw`\alpha_t=A_t\vee\mathbf{1}\{\mathrm{therapy}\wedge T_t\geq20000\},\quad\delta_g=1-(1-d_g)(1-\alpha_t k_g)`,
  'therapy_generation記錄首次啟動的當前generation。transition t到t+1使用剛更新的治療控制alpha。門檻以下若先前已治療仍繼續，因為therapy_generation不會重設None。自然存活與治療存活相乘，不能直接把兩種死亡率相加。基準敏感死亡delta=1-.7*.1=.93，抗藥delta=1-.7*.96=.328。');
 {
  const s=slide('狀態機：治療的歷史不能只看當前總量','相同負荷可能代表未治療，也可能已治療後下降',
   'A是當前門檻更新前的持續治療狀態，alpha為本次transition的治療狀態。完整狀態必須有A，因為下降到門檻以下後仍保持治療。停止模式C屬於程式控制，而非所有條件都是生物吸收。', 'modeling_cancer.py:214-244');
  const a=node(s,'未治療\nA = 0',110,245,300,110),b=node(s,'持續治療\nA = 1',850,245,300,110);connect(s,a,b);
  text(s,'therapy 開啟且 T ≥ 20,000',455,230,350,80,25, '#236B8E');
  text(s,'低於門檻時維持 A=0',110,400,350,70,26);text(s,'下降後仍維持 A=1',850,400,350,70,26);
  math(s,String.raw`S_t=(\mathbf{N}_t,A_t),\qquad A_{t+1}=\alpha_t`,70,505,1140,75);
  foot(s,'未停止時才執行轉移。停止條件與治療狀態分別控制觀察與生長。');
 }
 lesson('父代抽樣：兩次 Binomial 產生整批子代','一個基因型有 n 個父細胞，批次處理而不逐細胞迴圈',[[276,277],[285,291]],
  ['survivors = 存活父細胞 U_g\n以 n_g 和 1−δ_g 抽樣','triples = 三子代事件 V_g\n只在 survivors 中抽樣','n = 2U_g + V_g\n父代被新子代取代'],
  String.raw`U_g\sim\mathrm{Bin}(n_g,1-\delta_g),\quad V_g\mid U_g\sim\mathrm{Bin}(U_g,q_g),\quad B_g=2U_g+V_g`,
  '條件於state，父細胞死亡和分裂獨立同分布。第一次Binomial抽存活數，第二次在存活者中抽三子代數，故第二次試驗數是survivors而非count。每個存活者基準產生2個子代，三子代者再多1。U=0時V=0，沒有父代存活便不留下後代。每類使用批次Binomial與逐細胞Bernoulli在分布上等價。');
 {
  const s=slide('數值例子：基準細胞平均產生 1.428 個子代','先把存活和分裂機率結合，再理解生長或萎縮',
   '基準未治療delta=.3、q=.04。單個父細胞有0個子代的機率.3，2個子代.7*.96=.672，3個子代.7*.04=.028。均值0+2*.672+3*.028=1.428。治療下敏感均值.07*2.04=.1428，抗藥.672*2.04=1.37088，因而有選擇差異。期望不是一次試驗的必然細胞數。');
  table(s,[['子代 B','0','2','3'],['機率','0.300','0.672','0.028']],110,218,1060,155,null,31);
  math(s,String.raw`m_g(s)=\mathrm{E}[B_g^{(1)}\mid s]=(1-\delta_g(s))(2+q_g)`,70,407,1140,80);
  table(s,[['基準細胞','未治療','治療下敏感','治療下抗藥'],['每父細胞期望','1.428','0.1428','1.37088']],70,520,1140,90,null,25);
 }
 lesson('mutate_offspring：逐位元拆分子代計數','每輪只替尚未突變的位元抽樣，允許同代多位元突變',[[169,180]],
  ['mutants ~ Binomial(count, μ)\n一次抽出該類的突變細胞數','genotype | bit 設定新位元\nget(...,0) 合併相同基因型','依序遍歷全部 22 個位元\n本輪新類型下輪可再突變'],
  String.raw`X_{g,j}\sim\mathrm{Bin}(N_g,\mu),\qquad g\longmapsto g\cup\{j\}`,
  'population在函數入口先複製，updated在每一個gene開始時再複製。從舊population讀，向updated寫，避免本輪新基因型再處理同一個位元。當進入下一個gene時已接收前一輪更新，因此可在同一代新增兩個或更多不同位元。已為1的位元跳過，不允許回復突變。算法生成與每個子代每個未突變位元獨立Bernoulli相同的聯合分布。');
 lesson('突變守恆：重新分類細胞，不創造額外細胞','舊類型移出 X 個，新類型加入同樣的 X 個',[[181,187]],
  ['例：100 個 g=0，抽到 2 個突變\n結果為 {0:98, 8:2}','刪除 count=0 的項目\n維持稀疏字典與正計數','守恆只適用於突變階段\n存活與分裂會改變總數'],
  String.raw`\sum_hN_h^{\rm after\ mutation}=\sum_gB_g,\qquad g\subseteq h`,
  '例子抽樣數2是教學示意，不是由mu=.001必然得到。當gene=3遮罩為8，2個細胞由key0移至key8，總數100不變。dictionary.get容許新key存在或不存在。最後移除零數目。測試mu=0時完全不變，mu=1時全部移到滿位元基因型。');
 {
  const s=slide('逐位元抽樣對應精確突變矩陣 Q','Q 的列是父代基因型，欄是子代基因型',
   '把基因型表示為已突變位元集合g與h。只有g是h子集時可達。新增位元有|h|-|g|個，各以mu突變；保持0的L-|h|個位元各以1-mu不突變。已為1的位元不需再抽。當g不是h子集時Qgh=0。對任意g每列和1。給定Bg，Ygh聯合為Multinomial，這是逐位元拆分的等價數學表達，程式沒有真的建立Q。', 'modeling_cancer.py:147-169');
  math(s,String.raw`Q_{gh}=\mu^{|h|-|g|}(1-\mu)^{22-|h|}\quad(g\subseteq h),\qquad Q_{gh}=0\quad(g\nsubseteq h)`,70,222,1140,80);
  math(s,String.raw`(Y_{gh})_h\mid B_g\sim\mathrm{Multinomial}\bigl(B_g,(Q_{gh})_h\bigr),\qquad N_{h,t+1}=\sum_gY_{gh}`,70,353,1140,90);
  line(s,'程式逐位元抽樣，數學可用多項分布描述相同分布',485,75,1130,28);
  foot(s,'L=22。Q 每列總和為 1；它描述單個子代的突變，不是整個腫瘤轉移。');
 }
 {
  const s=slide('期望矩陣 M：繁殖能力乘上突變機率','M 記錄期望子代數，列總和為 m_g',
   '簡化只保留敏感S與抗藥R且關閉其他優勢。Q=[[.999,.001],[0,1]]，治療下繁殖均值mS=.1428、mR=1.37088，所以M=diag(m)Q。第一列和.1428，第二列和1.37088，並不是機率矩陣。這個兩類簡化例是教學縮減，不能當成原模型所有基因型的兩類精確閉合。', 'modeling_cancer.py:238-253');
  math(s,String.raw`M_{gh}(s)=m_g(s)Q_{gh},\qquad\mathrm{E}[\mathbf{N}_{t+1}\mid S_t=s]=\mathbf{n}M(s)`,70,212,1140,80);
  table(s,[['治療下 M','敏感 S','抗藥 R'],['敏感 S','0.1426572','0.0001428'],['抗藥 R','0','1.37088']],140,334,1000,190,null,29);
  math(s,String.raw`(100,5)M_{\rm on}=(14.26572,\ 6.86868)`,70,540,1140,65);
 }
 {
  const s=slide('完整 Markov 狀態：細胞組成加治療開關','給定 Sₜ，下一代分布不再需要更早的族群歷史',
   'N是所有基因型計數向量，A保存是否已開始治療，Config固定。IL11可由N計算，無需額外加入狀態。只用總細胞數T通常不具Markov性，因為相同T但抗藥比例或driver比例不同會有不同轉移。這裡的Markov鏈是理想隨機模型，實作中seed讓偽亂數路徑可重現。程式停止控制與復發歷史判讀需另外擴充狀態。');
  math(s,String.raw`S_t=(\mathbf{N}_t,A_t),\quad P_{ss'}=\Pr(S_{t+1}=s'\mid S_t=s)`,70,235,1140,80);
  math(s,String.raw`\Pr(S_{t+1}\mid S_t,S_{t-1},\ldots)=\Pr(S_{t+1}\mid S_t)`,70,365,1140,80);
  line(s,'相同 1,000 個細胞，抗藥組成不同，下一代分布不同',495,75,1130,28);
  foot(s,'固定 Config 時，完整 P 可固定。程式直接抽樣，不枚舉所有完整狀態。');
 }
 {
  const s=slide('完整 P 的數學表示：單個父細胞的生成函數','死亡、二子代、三子代與突變，合併成一個轉移規則',
   'z_h是子基因型h的形式變數。u_g(z)=sumQgh zh是一個子代的生成函數。單個父細胞死亡時為常數delta，產生2或3個獨立子代時分別為u²和u³，依對應概率加權得到Fg。這同時捕捉子代數與突變種類的聯合分布。程式以Binomial分批實作這個規則。');
  math(s,String.raw`u_g(\mathbf{z})=\sum_hQ_{gh}z_h`,70,217,1140,80);
  math(s,String.raw`F_g(\mathbf{z}\mid s)=\delta_g+(1-\delta_g)\left[(1-q_g)u_g(\mathbf{z})^2+q_gu_g(\mathbf{z})^3\right]`,70,352,1140,90);
  line(s,'δ_g：死亡　　u²：兩個子代　　u³：三個子代',505,75,1130,28);
  foot(s,'公式把 simulate() 與 mutate_offspring() 的抽樣規則一起寫下來。');
 }
 {
  const s=slide('完整 P 的數學表示：全族群生成函數','父細胞條件獨立，乘積的指定係數就是轉移機率',
   '給定s，ng個同基因型父細胞獨立，故該類生成函數Fg^ng。各類再相乘得到H。取z向量對應nprime的多項式係數得到下一代完整計數概率；另外要求下一代治療開關aprime等於本次更新後alpha。這個精確核定義於生物轉移的未停止狀態。停止控制可擴充狀態(N,A,t,C)。');
  math(s,String.raw`H_s(\mathbf{z})=\prod_g F_g(\mathbf{z}\mid s)^{n_g}`,70,220,1140,80);
  math(s,String.raw`P_{ss'}=\mathbf{1}\{a'=\alpha(s)\}\,[\mathbf{z}^{\mathbf{n}'}]H_s(\mathbf{z})`,70,355,1140,80);
  math(s,String.raw`[\mathbf{z}^{\mathbf{n}'}]\ :\ \mathrm{coefficient\ of}\ \prod_h z_h^{n'_h}`,70,487,1140,70);
  foot(s,'這是解析定義。直接儲存如此巨大的 P，成本遠高於批次抽樣。');
 }
 {
  const s=slide('P 的數值例：把狀態縮小後就能列出矩陣','停用突變與治療，N≥4 的首次進入以 † 記帳',
   '狀態順序0,1,2,3,dagger。基準單父子代0,2,3機率.3,.672,.028。例如2個父細胞都死亡=.09，一死一產2=.4032，一死一產3=.0168，其餘進入overflow概率.49。3個父細胞對應.027,.18144,.00756,.784。dagger是首次N>=4後不再觀察的人工吸收，不是環境飽和或完整模型中的生物吸收。每列總和1。');
  table(s,[['列／欄','0','1','2','3','†'],...P.map((r,i)=>[['0','1','2','3','†'][i],...r.map(v=>v===0||v===1?String(v):v.toFixed(5))])],95,215,1090,325,null,28);
  math(s,String.raw`P_{ij}=\Pr(S_{t+1}=j\mid S_t=i),\qquad\sum_jP_{ij}=1`,70,555,1140,60);
 }
 {
  const s=slide('為何 π 乘在 P 左側：列向量的加權平均','每個下一步狀態，累加所有來源狀態的貢獻',
   '本簡報採橫向row probability vector。pi_t(i)是現在在狀態i的概率，Pij是由i走到j的條件概率，對所有i加總得到pi_next(j)。因此維度(1×K)(K×K)=(1×K)。若改column convention可寫p_next=Ptranspose p，不是左右哪個固定正確，而是索引與向量方向必須一致。pi0=[0,1,0,0,0]表示初始確定在一個細胞，pi1=P第二列。');
  math(s,String.raw`\pi_{t+1}(j)=\sum_i\pi_t(i)P_{ij},\qquad\boldsymbol{\pi}_{t+1}=\boldsymbol{\pi}_tP`,70,220,1140,90);
  math(s,String.raw`(1\times K)(K\times K)=(1\times K)`,70,353,1140,80);
  math(s,String.raw`\boldsymbol{\pi}_0=(0,1,0,0,0),\quad\boldsymbol{\pi}_1=(0.30,0,0.672,0.028,0)`,70,482,1140,85);
  foot(s,'π 的項目是「整個族群處於各狀態的機率」，不是各基因型細胞比例。');
 }
 {
  const s=slide('矩陣連乘描述所有路徑，不是一條模擬路徑','P² 加總經過所有中間狀態的兩步轉移機率',
   '從一個細胞兩代後滅絕概率=.3*1+.672*.09+.028*.027=.361236。每一項是經過0、2、3個細胞中間狀態的路徑概率。完整pi2由pi0P²計算且總和1。一條simulate路徑只抽中其中一些狀態；P連乘描述全部可能路徑的分布。路徑probability初始pi0(s0)乘沿途P，不將實際抽中的計數向量直接乘P。');
  math(s,String.raw`(P^2)_{ij}=\sum_kP_{ik}P_{kj},\qquad\boldsymbol{\pi}_t=\boldsymbol{\pi}_0P^t`,70,217,1140,90);
  math(s,String.raw`\Pr(N_2=0)=0.30+0.672(0.09)+0.028(0.027)=0.361236`,70,360,1140,85);
  math(s,String.raw`\Pr(s_0,\ldots,s_T)=\pi_0(s_0)\prod_{t=0}^{T-1}P_{s_t s_{t+1}}`,70,502,1140,80);
  foot(s,'玩具例子的數字不能直接取代完整腫瘤的 P。單條路徑使用沿途機率乘積。');
 }
 {
  const s=slide('P、Q 與 M：三種矩陣回答不同問題','索引不同，因此不能把它們互相替代',
   'P兩個索引是完整族群狀態，元素機率，推進pi。Q索引是單細胞父與子基因型，元素突變概率。M索引也是基因型，但元素期望子代數，列和m，不保證1。固定環境時E[Nt]=N0M^t；當IL11與治療由隨機state決定時E[Nnext]=E[Nt M(St)]，一般不等於E[Nt] M(E[St])。完整P可固定，因不同state的環境寫在不同列。');
  table(s,[['矩陣','索引','元素','列總和'],['P','完整族群狀態','下一步機率','1'],['Q','父／子基因型','子代突變機率','1'],['M','父／子基因型','期望子代數','m_g']],70,215,1140,245,[190,365,385,200],27);
  math(s,String.raw`\mathrm{E}[\mathbf{N}_{t+1}]=\mathrm{E}[\mathbf{N}_tM(S_t)],\quad\mathrm{E}[\mathbf{N}_t]=\mathbf{N}_0M^t\ \mathrm{only\ if}\ M\ \mathrm{is\ fixed}`,70,505,1140,85);
  foot(s,'IL11 與治療使 M 隨當前狀態改變。完整 Markov 核 P 仍可保持固定。');
 }
 lesson('summarize：把基因型字典變成可畫圖的統計','統計只讀當前族群，與下一代抽樣分開',[[192,196]],
  ['total：總細胞數\nresistant：抗藥位元為 1 的數量','richness：不同基因型數\nShannon：考慮比例與均勻度','滅絕時比例避免除以 0\n細胞數與比例的定義分開'],
  String.raw`p_g=\frac{N_g}{T},\quad H=-\sum_gp_g\log p_g,\quad D=1-\sum_gp_g^2`,
  'Shannon用自然對數。richness=len(population)只看種類，不看大小比例。例{a:50,b:50}的H=ln2≈.693，而{a:99,b:1}的H≈.056，比較相同richness但均勻度不同。exp(H)稱有效基因型數。零族群Shannon=0且effective_genotypes另外定為0，避免exp(0)=1誤表示仍有一種細胞。基因型多樣性不是家系多樣性。');
 lesson('therapy_outcome：從 history 判讀復發事件','先低於治療起始負荷，再回到該負荷才算復發',[[301,308],[310,314]],
  ['baseline 是開始治療時總量\n不是固定寫死的 20,000','declined 保存「曾下降」\nrelapse 是第一次回升代數','未啟動、未觀察到復發與滅絕\n保留不同 outcome 標記'],
  String.raw`T_*=T_{t_{\rm therapy}},\quad t_{\rm relapse}=\inf\{t>t_{\rm therapy}:T_t\geq T_*,\ \exists u<t:T_u<T_*\}`,
  'therapy_generation對應history的generation索引，因從0開始每代記錄。baseline可overshoot20,000。after排除啟動代。declined一旦true保持true，防止單純持續增長被叫復發。完整分類依序relapse、extinction、no_relapse_observed、no_initial_response；未開始治療另為therapy_not_started。若要在Markov狀態內實時識別這個歷史事件，需加入baseline與declined，不只是N,A。');
 lesson('main：用獨立試驗建立可重現的實驗流程','同一 seed 重現同一路徑，不同 replicate 使用不同 seed',[[781,784],[793,794],[797,797]],
  ['治療組 seed=42+100000+rep\nrep=0,…,29，共 30 次','Config、seed、套件版本\n與 outcome 一起寫入 metadata','CSV 保留逐代和逐次資料\n簡報圖可回溯到原始數字'],
  String.raw`S_{t+1}^{(r)}\mid S_t^{(r)}=s\sim P(s,\cdot),\quad r=1,\ldots,30`,
  '相同參數、seed與相同執行環境通常重現該偽亂數路徑。不同trial不共用同一Generator，避免意外耦合。main另外執行one/twohit與IL11組，比較時不施加治療。seed是重現標記，不是模型狀態或生物參數。套件版本存在metadata，有助查詢環境變化。');
 lesson('多次試驗統計：滅絕零值與停止缺失分開','後期曲線必須同時知道有效試驗數 n',[[694,694],[698,701],[713,714]],
  ['資料矩陣先填 NaN\n只寫入實際觀察或滅絕延續','計算上限之後維持 NaN\n滅絕之後細胞數延續 0','比例只對仍活且有觀察的腫瘤\n每代保留各指標有效 n'],
  String.raw`\widehat{Q}_{p,t}=\mathrm{quantile}_p\{X_t^{(r)}:\ X_t^{(r)}\ \mathrm{observed}\},\quad p\in\{0.1,0.5,0.9\}`,
  'plot_therapy_ensemble依絕對generation對齊，不是對齊治療起始。停止於celllimit的trial後續是缺失，不能填0或延續最後值。extinction後細胞數為0；抗藥比例及基因頻率無分母，因此缺失。percentiles是trial間變異，不是confidenceinterval。多種曲線的中位數不必相加。有效n下降時需小心條件選擇。');
 lesson('等大小取樣：避免跨門檻大小不同的混淆','從首次達標的完整族群無放回抽取 5,000 個細胞',[[425,433]],
  ['不同死亡率達標代數不同\n最後族群可能大於 5,000','multivariate_hypergeometric\n保留固定樣本總數且不放回','比較等大小樣本的 H\n仍須注意達標時間與選擇'],
  String.raw`\mathbf{X}\sim\mathrm{MultiHypergeom}(\mathbf{N},5000),\qquad\sum_gX_g=5000`,
  'sample_at_size從最後達標族群取sample，不是連續時間精確N=5000時的整個腫瘤。超幾何分布是無放回抽樣，與multinomial的有放回近似不同。未達標trial回None並記錄。高死亡可能較慢而累積更多mutation機會，故高H不能單獨歸因於death造成多樣性。下一步固定時間與固定大小分別比較。');
 {
  const s=slide('test_modeling_cancer.py：測試生物語意與不變量','測試回答「程式是否實作了預定的模型」',
   '現有五個unittest方法涵蓋IL11共享、ensemble缺失與滅絕、同一TSG兩擊、突變守恆與獨立性、seed重現及參數驗證。mu=.2只有一個TSG兩位元可突變時，0/1/2擊概率.64/.32/.04，10萬樣本比例容許.006誤差。測試能抓算法錯誤，但不能證明教學參數符合生物數據。', 'test_modeling_cancer.py: SuppressorTests');
  table(s,[['檢查','期望','防止的錯誤'],['TSG 語意','不同基因各一擊不算兩擊','把 one-hit 和 two-hit 混淆'],['突變守恆','前後總細胞數相同','分類時重複或遺失細胞'],['獨立等位基因','0/1/2 擊約 0.64/0.32/0.04','同代多突變處理錯誤'],['IL11、截尾、seed','共享效益、缺失、重現','環境時序與統計語意錯誤']],70,215,1140,340,[230,455,455],25);
  math(s,String.raw`K\sim\mathrm{Bin}(2,0.2):\quad\Pr(K=0,1,2)=(0.64,0.32,0.04)`,70,565,1140,55);
 }
 {
  const s=slide('實驗設計：每張結果圖對應不同控制條件','程式的參數與取樣方式決定可以比較什麼',
   '主治療30次用twohit且IL11off。onehit/twohit18代無治療保持相同loci與oncogene，只改TSG門檻且不同seed組。IL11on/off18代無治療相同loci及同一replicate seed，族群分歧後不代表逐細胞配對。deathrate.2/.35/.45首次達5000後固定樣本，達標代數不同。這些比較估計的是指定教學參數情況的行為。', 'modeling_cancer.py: run_comparisons, run_il11_comparison, main');
  table(s,[['問題','改變條件','固定條件／比較方式'],['治療與復發','族群達門檻後啟動','two-hit、IL11 OFF、30 次'],['TSG 門檻','one-hit / two-hit','無治療、18 代、每組 30 次'],['IL11 共享','ON / OFF','無治療、同位點、每組 30 次'],['多樣性','d=0.20 / 0.35 / 0.45','每次取 5,000 細胞、每組 30 次']],70,215,1140,350,[235,395,510],25);
  foot(s,'下一頁開始以實際模擬結果檢查上述程式規則的表現。');
 }
}
