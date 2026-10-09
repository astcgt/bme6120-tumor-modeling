/** Rebuild the Traditional Chinese teaching deck with the bundled Artifact Tool.
 * Node runtime and Python may be overridden through RUNTIME_NODE_MODULES and
 * RUNTIME_PYTHON. No project Python or conda dependencies are changed.
 */
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import {pathToFileURL} from 'node:url';
import {createRequire} from 'node:module';

const workspace = path.resolve(import.meta.dirname, '..');
const runtime = path.join(os.homedir(), '.cache/codex-runtimes/codex-primary-runtime/dependencies');
const modules = process.env.RUNTIME_NODE_MODULES ?? path.join(runtime, 'node/node_modules');
process.env.RUNTIME_NODE_MODULES = modules;
const req = createRequire(path.join(modules, '__runtime__.cjs'));
const {Presentation, PresentationFile, FileBlob} = await import(pathToFileURL(req.resolve('@oai/artifact-tool')).href);
const skill = path.join(os.homedir(), '.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations');
const {finalizePresentation,applyPresentationChartFont} = await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const build = process.env.PRESENTATION_BUILD_DIR ?? path.join(workspace, '.presentation-build');
const out = path.join(workspace, 'outputs/presentations');
await fs.mkdir(build,{recursive:true});
await fs.mkdir(out,{recursive:true});

const INK='#173042', BLUE='#236B8E', RED='#A74F3D', GRAY='#657580', BG='#FCFBF8';
const FONT='Arial Unicode MS';
const deck=Presentation.create({slideSize:{width:1280,height:720}});
const records=[];
const nativeTables=[];
const chartOwners=[];
const results=JSON.parse(await fs.readFile(path.join(workspace,'presentations/data/results.json'),'utf8'));

function text(s,str,x,y,w,h,size=28,color=INK,bold=false,align='left'){
  const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  sh.text=str;
  sh.text.style={typeface:FONT,fontSize:size,color,bold,alignment:align,verticalAlignment:'top',autoFit:'none',wrap:'square',insets:{left:0,right:0,top:0,bottom:0}};
  return sh;
}
function slide(title,subtitle,notes,source='modeling_cancer.py: Config, simulate, mutate_offspring'){
  const s=deck.slides.add();s.background.fill=BG;
  const n=records.length+1;
  text(s,title,64,44,1152,72,43,INK,true);
  if(subtitle)text(s,subtitle,66,125,1148,66,26,GRAY);
  text(s,String(n).padStart(2,'0'),1170,683,48,26,18,GRAY,false,'right');
  s.speakerNotes.textFrame.setText(notes+'\n\n來源：本專案 '+source+'。參數為教學假設，未用臨床資料校準。');
  records.push({number:n,title,subtitle,notes,source});
  return s;
}
function eq(s,str,y=235,size=38,x=70,w=1140,h=75){return text(s,str,x,y,w,h,size,BLUE);}
function line(s,str,y,x=70,w=1140,size=28){return text(s,str,x,y,w,74,size);}
function foot(s,str){text(s,str,68,619,1100,60,22,GRAY);}
function node(s,str,x,y,w=260,h=90,color=BLUE){
  const sh=s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:'#F0F5F7',line:{fill:color,width:1.5}});
  sh.text=str;sh.text.style={typeface:FONT,fontSize:27,color:INK,alignment:'center',verticalAlignment:'middle',autoFit:'none',insets:{left:9,right:9,top:5,bottom:5}};
  return sh;
}
function connect(s,a,b,from='right',to='left',kind='straight'){
  return s.shapes.connect(a,b,{kind,fromSide:from,toSide:to,line:{fill:GRAY,width:2},tail:{type:'triangle',length:'med',width:'med'}});
}
function table(s,values,x,y,w,h,widths=null,size=27){
  const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,values:values.map(r=>r.map(String)),...(widths?{columnWidths:widths}:{})});
  nativeTables.push(records.length);
  t.styleOptions={headerRow:false,bandedRows:false};
  t.borders.assign({fill:'#D9E1E5',width:.6,style:'solid'});
  for(let i=0;i<values.length;i++)for(let j=0;j<values[0].length;j++){
    const c=t.getCell(i,j);c.fill=i===0?'#E8F0F4':BG;
    c.text.style={typeface:FONT,fontSize:size,color:INK,bold:i===0,alignment:'center',verticalAlignment:'middle',autoFit:'none',insets:{left:6,right:6,top:5,bottom:5}};
  }
  return t;
}
const P=[[1,0,0,0,0],[.3,0,.672,.028,0],[.09,0,.4032,.0168,.49],[.027,0,.18144,.00756,.784],[0,0,0,0,1]];
function vmul(v,m){return m[0].map((_,j)=>v.reduce((sum,z,i)=>sum+z*m[i][j],0));}
const pi0=[0,1,0,0,0],pi1=vmul(pi0,P),pi2=vmul(pi1,P);
for(const row of P)if(Math.abs(row.reduce((a,b)=>a+b,0)-1)>1e-12)throw Error('P row sum');
if(Math.abs(pi2[0]-.361236)>1e-12)throw Error('P squared');

function scientificChart(s,kind,categories,series,ytitle,{xTitle='世代',min=0,max=undefined,labels=false}={}){
 chartOwners.push(records.length);
 const ch=s.charts.add(kind,{
  position:{left:65,top:200,width:875,height:380},categories:categories.map(String),series,
  hasLegend:kind!=='bar',legend:{position:'bottom',textStyle:{typeface:FONT,fontSize:19,fill:INK}},
  chartFill:BG,chartLine:{fill:'none',width:0},plotAreaFill:BG,plotAreaLine:{fill:'none',width:0},
  xAxis:{title:{text:xTitle,textStyle:{typeface:FONT,fontSize:22,fill:INK}},textStyle:{typeface:FONT,fontSize:19,fill:INK},line:{fill:GRAY,width:1}},
  yAxis:{title:{text:ytitle,textStyle:{typeface:FONT,fontSize:22,fill:INK}},min,...(max?{max}:{}),numberFormatCode:'#,##0',
    textStyle:{typeface:FONT,fontSize:19,fill:INK},majorGridlines:{fill:'#E0E5E8',width:.6},line:{fill:GRAY,width:1}},
  lineOptions:{smooth:false},scatterOptions:{style:'marker'},
  barOptions:{direction:'column',grouping:'clustered',gapWidth:80},
  dataLabels:{showValue:labels,position:'outEnd',textStyle:{typeface:FONT,fontSize:22,fill:INK}},
 });
 applyPresentationChartFont(ch,{fontFamily:FONT});
 return ch;
}
function curve(name,values,color,percentile=false){return {name,values,fill:color,
 line:{fill:color,width:percentile?1.3:3,style:percentile?'dashed':'solid'},marker:{symbol:'none'}};}
function side(s,stat,body){text(s,stat,970,235,240,150,38,BLUE,true);text(s,body,970,408,240,180,24);}

// 01
{
 const s=slide('腫瘤演化的 Markov 模型','隨機細胞繁殖、狀態轉移與矩陣連乘',
 '這份簡報先用一個細胞的繁殖理解隨機模型，再引入完整狀態和機率矩陣，接著呈現實際模擬結果與後續實驗。最後四頁為數學附錄，可依聽眾背景補充。模型採同步世代，時間沒有換算成天數。注意始終區分單次路徑、機率分布及期望細胞數。');
 text(s,'Sₜ = (Nₜ, Aₜ)',78,252,1100,100,66,BLUE);
 text(s,'πₜ = π₀ Pᵗ',78,382,1100,100,66,BLUE);
 text(s,'BME 6120　繁體中文教學簡報',80,550,1080,54,29,GRAY);
}
// 02
{
 const s=slide('模型觀察的是整個細胞族群','初始 100 個未突變細胞，逐代追蹤族群組成',
 '程式不追蹤每個細胞的家系，而是把相同基因型的細胞合併計數。初始族群由一百個零突變細胞構成。每代所有父細胞由零、二或三個子細胞取代。詢問觀眾：若某種突變使細胞留下更多後代，族群比例會如何改變？治療會進一步改變選擇條件。');
 const a=node(s,'當代族群\n基因型與細胞數',70,265,300,120);
 const b=node(s,'死亡與繁殖\n子代再發生突變',490,265,300,120);
 const c=node(s,'下一代族群\n新數量與新組成',910,265,300,120);
 connect(s,a,b);connect(s,b,c);
 line(s,'記錄總負荷、抗藥比例與基因型多樣性',465);
 foot(s,'同一基因型可來自不同祖先，基因型計數不等於家系重建。');
}
// 03
{
 const s=slide('一個父細胞的三種結果','死亡機率 d = 0.30，存活後三子細胞機率 q = 0.04',
 '先死亡或存活，再讓存活者產生二或三個子細胞。產生三子細胞的無條件機率是零點七乘零點零四，等於零點零二八。請觀眾檢查三個機率加總為一。父細胞不會額外留在下一代。期望一點四二八不是每次必定得到的數量。', 'modeling_cancer.py:49-60, 245-249');
 const a=node(s,'1 個父細胞',76,330,260,100);
 for(const [y,str,p] of [[210,'0 個子細胞','d = 0.300'],[355,'2 個子細胞','(1−d)(1−q) = 0.672'],[500,'3 個子細胞','(1−d)q = 0.028']]){
   const b=node(s,`${str}\n${p}`,650,y,540,110);connect(s,a,b);
 }
 foot(s,'E[X] = 0×0.300 + 2×0.672 + 3×0.028 = 1.428');
}
// 04
{
 const s=slide('突變改變細胞的繁殖能力','Oncogene 與 TSG 優勢可相加，每一類最多計算一次',
 '預設 oncogene 任一個突變即達門檻。TSG 預設要求同一個基因的兩個等位基因失活，兩個不同 TSG 各一擊不符合 two-hit。表中的 d 不包含治療，IL11 關閉。平均子代數是單個父細胞的繁殖能力指標。m 大於一表示期望擴張，但小族群仍可能隨機滅絕。', 'modeling_cancer.py:105-132');
 table(s,[['表型','自然死亡 d','三子細胞 q','E[X]'],['未突變','0.30','0.04','1.4280'],['Oncogene 達標','0.22','0.16','1.6848'],['Oncogene 與 TSG 達標','0.14','0.28','1.9608']],70,225,1140,280,[440,230,230,240]);
 eq(s,'m = E[X] = (1−d)(2+q)',535,37);
 foot(s,'TSG two-hit 要求同一基因雙等位失活。所有數值均為模型假設。');
}
// 05
{
 const s=slide('完整狀態包含族群組成與治療開關','知道現在的完整狀態，就能決定下一代的機率分布',
 'N 粗體在數學表示中是一個向量，每一項是某種基因型的細胞數，程式以 population 字典儲存。A 是本代門檻檢查之前是否已經啟動治療。過去的影響透過累積突變和治療開關進入現在。當新門檻觸發治療，這次轉移使用更新後的 alpha。', 'modeling_cancer.py:196-244');
 eq(s,'Sₜ = (Nₜ, Aₜ)',230,54);
 line(s,'Nₜ = (N₀,ₜ, N₁,ₜ, …)　各基因型的細胞數',350);
 line(s,'Aₜ ∈ {0,1}　治療尚未啟動／已啟動',435);
 line(s,'例：Sₜ = ((100 個敏感，5 個抗藥), 1)',525,70,1140,31);
 foot(s,'此例只保留兩種類型。完整模型使用 22 個位元表示 20 個基因座。');
}
// 06
{
 const s=slide('只知道總細胞數，資訊仍不足','同樣 1,000 個細胞，治療後的結果可以相差很大',
 '此頁用簡化的兩種類型說明為何總負荷本身通常不具有 Markov 性質。假設所有細胞自然 d 為零點三、q 為零點零四，沒有 driver、TSG 或 IL11 優勢，且治療正在進行。敏感 m 等於零點一四二八，抗藥 m 等於一點三七零八八。因此 A 的期望總量為一百五十五點零八零八，B 為一千三百五十八點五九九二。即使敏感與抗藥數量相同，其他突變仍可影響完整模型。');
 table(s,[['族群','敏感','抗藥','治療後期望總數'],['A','990','10','155.0808'],['B','10','990','1358.5992']],70,240,1140,230,[230,220,220,470]);
 eq(s,'E[Tₜ₊₁] = Nₛ×0.1428 + Nᵣ×1.37088',515,35);
 foot(s,'簡化計算：固定基準生長參數，停用其他生長優勢與 IL11。');
}
// 07
{
 const s=slide('程式狀態機與治療門檻','每代先記錄並檢查停止條件，再更新治療與推進一代',
 '這是程式控制模式，而不是完整族群 Markov 狀態。從初始未治療開始，每代先檢查停止。只有能繼續推進時才檢查治療門檻。未治療且首次達到兩萬時，下一次轉移立即施加治療。治療開始後即使負荷下降也不取消。兩個活動模式都能因停止條件結束。停止檢查的優先順序是滅絕、取樣目標、細胞上限、最後代數。復發屬於事後軌跡分類。', 'modeling_cancer.py:213-253, therapy_outcome');
 const a=node(s,'未治療\nA = 0',125,260,320,120);
 const b=node(s,'持續治療\nA = 1',820,260,320,120);
 connect(s,a,b);text(s,'治療功能開啟且 T ≥ 20,000',475,245,330,74,25,BLUE);
 const c=node(s,'停止\n滅絕／目標／上限／代數',390,480,500,100,RED);
 connect(s,a,c,'bottom','left','elbow');connect(s,b,c,'bottom','right','elbow');
 text(s,'未達條件時，繁殖與突變後\n仍留在未治療模式',80,390,405,74,24);
 text(s,'未達停止條件時持續治療\n負荷下降也不關閉',815,390,405,74,24);
 foot(s,'停止優先序：滅絕 > 取樣目標 > 計算上限 > 最後觀察代。');
}
// 08
{
 const s=slide('一次轉移的抽樣方程','對每個基因型 g，先抽存活者，再抽三子細胞事件',
 'delta 是包含治療的有效死亡機率，q 是存活後產生三子細胞的機率。給定當前完整狀態，不同父細胞的抽樣獨立。先抽 U，再在 U 條件下抽 V，因此兩者不是無條件獨立。子代數 B 等於二 U 加 V。程式最後再讓子代突變，新的生存優勢從下一次轉移開始生效。', 'modeling_cancer.py:237-253');
 eq(s,'U_g ∼ Binomial(n_g, 1−δ_g)',225,39,70,700);
 text(s,'存活父細胞數',790,238,400,55,27,GRAY);
 eq(s,'V_g | U_g ∼ Binomial(U_g, q_g)',340,39,70,700);
 text(s,'其中產生三個子細胞者',790,353,410,55,27,GRAY);
 eq(s,'B_g = 2U_g + V_g',455,43);
 line(s,'子代繼承父代突變，接著獨立新增突變',555,70,1140,28);
 foot(s,'δ 包含自然死亡與治療殺傷，q 是存活條件下的機率。');
}
// 09
{
 const s=slide('P 是完整族群狀態的機率表','橫列 row 表示現在，直欄 column 表示下一步',
 '此頁矩陣條目是符號，而非數值估計。每個 s 索引代表一整個族群及治療開關。P 的 row 是從該狀態出發的分布，所有可能下一步加總為一。請觀眾指出 P12 對應從 s1 到 s2，而不是從基因1到基因2。完整狀態集合極大，未加停止時為可數無窮。');
 eq(s,'Pᵢⱼ = Pr(Sₜ₊₁ = sⱼ | Sₜ = sᵢ)',225,39);
 table(s,[['現在／下一步','s₁','s₂','s₃','…'],['s₁','P₁₁','P₁₂','P₁₃','…'],['s₂','P₂₁','P₂₂','P₂₃','…'],['s₃','P₃₁','P₃₂','P₃₃','…']],160,330,960,225,[300,165,165,165,165]);
 foot(s,'Pᵢⱼ ≥ 0，且對每個出發狀態 i，∑ⱼ Pᵢⱼ = 1。');
}
// 10
{
 const s=slide('可以精確計算的五狀態範例','單一細胞類型，停用突變與治療，d = 0.30，q = 0.04',
 '此頁是一個明確定義的簡化分枝過程。狀態是零、一、二、三和首次達到四或更多之後的停止記帳狀態 dagger。從 i 個父細胞產生 j 個子細胞的機率可由生成函數 f 的 i 次方取係數。零是生物滅絕。dagger 的自迴圈是人為停止，不能把不受限制的所有大族群合併後仍假定 Markov。第一個狀態一代後不會仍只有一個細胞，因為每個父代只產生零、二或三個子代。', 'plot_markov_model.py: finite_count_matrix');
 table(s,[['現在／下一步','0','1','2','3','†'],...P.map((r,i)=>[i===4?'†':String(i),...r.map(v=>String(v))])],75,220,1130,310,[260,174,174,174,174,174],25);
 line(s,'例：P₂,₂ = 2×0.30×0.672 = 0.4032',553,80,1120,30);
 foot(s,'† 表示首次進入 N ≥ 4 後停止。每一列機率總和為 1。');
}
// 11
{
 const s=slide('π 記錄各狀態的機率','細胞數向量 N 與機率向量 π，回答不同問題',
 '此頁避免最常見的混淆。N 是實際細胞數，而 pi 是所有完整族群狀態的分布。五狀態範例從一個細胞開始，初始狀態確定，所以只有對應狀態一的項目為一。把 pi0 乘 P 之後得到一代的分布。一次隨機模擬只會抽到零、二或三個其中一個結果。');
 eq(s,'狀態順序：(0, 1, 2, 3, †)',225,35);
 eq(s,'π₀ = (0, 1, 0, 0, 0)',335,48);
 line(s,'目前確定有 1 個細胞，所以該狀態的機率是 1',425);
 eq(s,'π₁ = (0.300, 0, 0.672, 0.028, 0)',505,39);
 foot(s,'π 的每一項都是機率，其總和為 1。單次抽樣的細胞數仍是整數。');
}
// 12
{
 const s=slide('πP 是對所有出發狀態做加權','下一步到達 j 的機率，等於所有來源貢獻的總和',
 'pi 位於左邊是因為選擇橫向 row vector。pi 是一乘 K，P 是 K乘K，相乘得到一乘K。一般公式由全機率公式得到。若採直向 column vector，須使用 P 的轉置左乘。沒有任何生物學原因要求乘右邊，這是索引方向與向量形狀的約定。請先讀乘法的維度，再讀概率意義。');
 eq(s,'πₜ₊₁,ⱼ = ∑ᵢ πₜ,ᵢ Pᵢⱼ',220,43);
 text(s,'現在在 i 的機率',110,340,440,65,31);
 text(s,'從 i 到 j 的機率',675,340,470,65,31);
 eq(s,'πₜ₊₁ = πₜ P',445,48);
 line(s,'(1×K) × (K×K) = (1×K)',532,75,1110,31);
 foot(s,'若使用直向向量 pₜ，則 pₜ₊₁ = Pᵀpₜ。兩種寫法表示同一個過程。');
}
// 13
{
 const s=slide('兩代後滅絕的機率','P² 自動加總所有可能的中間路徑',
 '沿用五狀態範例且從一個細胞開始。第一代可能為零、二或三。第二代滅絕機率為零點三加零點六七二乘零點零九加零點零二八乘零點零二七，等於零點三六一二三六。這裡三條路徑互斥，所以可相加。單條路徑概率乘積和多條路徑概率和是兩種操作。');
 table(s,[['中間路徑','路徑機率'],['1 → 0 → 0','0.300 × 1 = 0.300000'],['1 → 2 → 0','0.672 × 0.090 = 0.060480'],['1 → 3 → 0','0.028 × 0.027 = 0.000756']],85,220,1110,280,[465,645],29);
 eq(s,'Pr(N₂ = 0) = 0.361236',540,43);
 foot(s,'一般形式：(P²)ᵢⱼ = ∑ₖ PᵢₖPₖⱼ。');
}
// 14
{
 const s=slide('多代演化可以用矩陣連乘','固定參數與完整狀態下，使用同一個轉移矩陣',
 '在完整狀態包含 N 和 A 時，治療與IL11狀態依賴已經進入每一列轉移規則，因此基本族群過程可時間齊次。pi0 P平方保留所有結果分布。若外部參數由時間決定才使用P0P1。固定最後觀察代數是觀察停止時間，可另納入時鐘以建立包含控制停止的鏈。細胞死亡和繁殖的隨機性並不使P本身必須每代改變。');
 eq(s,'πₜ = π₀ Pᵗ',215,58);
 table(s,[['世代 2 的狀態','0','1','2','3','†'],['π₂','0.361236','0','0.276031','0.011501','0.351232']],75,345,1130,145,[260,174,174,174,174,174],25);
 line(s,'若參數隨時間改變：πₜ = π₀ P₀P₁…Pₜ₋₁',530,70,1140,31);
 foot(s,'表中為六位小數顯示。完整模型的 P 極大，程式透過抽樣產生路徑。');
}
// 15
{
 const s=slide('Q 描述子細胞的基因型突變','只保留抗藥位元的簡化例子，μ = 0.001',
 'Q 的列代表父代基因型，欄代表子代基因型，條件是在子細胞已經產生之後。敏感子代以零點零零一機率獲得抗藥突變。抗藥不回復。死亡不在Q內，因為死亡由父細胞繁殖分布處理。完整模型Q共有2的22次方潛在基因型，且允許同一代多位元突變。');
 table(s,[['父代／子代','敏感 S','抗藥 R'],['敏感 S','0.999','0.001'],['抗藥 R','0','1']],210,240,860,240,[330,265,265],34);
 line(s,'每一列總和為 1，已突變位元保持突變',530,75,1140,31);
 foot(s,'Q 的索引是基因型。P 的索引是整個族群狀態。');
}
// 16
{
 const s=slide('M 描述每個父細胞的期望子代','繁殖能力 m 乘上子代突變機率 Q',
 'Mgh等於一個g型父細胞的期望子代數乘上子代為h的機率。因此M不是row stochastic，其列和是m。本頁與後兩頁均為只保留抗藥位元、其他優勢關閉的簡化計算。未治療敏感和抗藥m都為1.428，治療後分別0.1428和1.37088。新抗藥突變在子代產生後發生，敏感父細胞仍按敏感殺傷率。');
 eq(s,'M_g,h(s) = m_g(s) Q_g,h',215,44);
 eq(s,'m_g(s) = (1−δ_g(s))(2+q_g)',310,40);
 table(s,[['細胞類型','未治療 m','治療下 m'],['基準敏感','1.428','0.1428'],['基準抗藥','1.428','1.37088']],95,410,1090,178,[430,330,330],29);
 foot(s,'治療殺傷作用於自然存活者：敏感 k = 0.90，抗藥 k = 0.04。');
}
// 17
{
 const s=slide('治療下的期望矩陣與計算','目前有 100 個敏感細胞與 5 個抗藥細胞',
 '這裡以兩種類型順序S,R顯示Mon。第一列總和0.1428，第二列1.37088。下一代抗藥期望6.86868包含5個抗藥父細胞貢獻6.8544，加上100個敏感父細胞貢獻0.01428個新抗藥子代。期望值可以小於一個，表示大量重復試驗平均。一次試驗不能產生分數個細胞。');
 table(s,[['父代／子代','敏感 S','抗藥 R'],['敏感 S','0.1426572','0.0001428'],['抗藥 R','0','1.37088']],130,220,1020,215,[390,315,315],31);
 eq(s,'E[Nₜ₊₁ | Sₜ=s] = (100,5) M_on',470,36);
 eq(s,'= (14.26572, 6.86868)',550,40);
}
// 18
{
 const s=slide('P、Q 與 M 的索引和用途','同樣是矩陣，但各自記錄不同的量',
 '請讓觀眾說明三個矩陣各自一個元素的意義。P是完整族群狀態間的機率，Q是已產生子代的突變機率，M是單個父細胞的期望各類子代數。P可推算完整分布、滅絕概率等，M只給條件期望。Q對死亡沒有預測，因為先條件在子代已經出生。');
 table(s,[['矩陣','索引代表','元素的意義','每列總和'],['P','完整族群狀態','下一步狀態機率','1'],['Q','父／子基因型','子代突變機率','1'],['M','父／子基因型','期望子代數','m_g']],70,235,1140,285,[170,320,400,250],28);
 line(s,'機率分布：πₜ₊₁ = πₜP',552,75,530,30);
 text(s,'細胞數期望：E[N′|s] = nM(s)',620,552,580,68,30,BLUE);
}
// 19
{
 const s=slide('IL11 讓下一代依賴當前族群組成','生產者改善共享環境，非生產者也得到相同效益',
 '開啟IL11後位元6由中性passenger變為生產者位點。f是當前生產者比例，b為全體細胞共享的死亡率降低量。百分之二生產者給零點零三絕對降低，也就是三個百分點。不是相對下降百分之三。對N向量完整狀態，f和b是確定函數，所以不破壞Markov性質。若加入細胞因子殘留，要把當前濃度也放入狀態。', 'modeling_cancer.py:135-144');
 eq(s,'f = IL11⁺ 細胞數 / 總細胞數',215,35);
 eq(s,'b(N) = 0.06 f / (0.02+f)',310,43);
 table(s,[['IL11⁺ 比例 f','共享降低量 b','基準自然死亡率'],['0%','0','30%'],['2%','0.03','27%'],['10%','0.05','25%']],110,415,1060,180,[350,350,360],28);
 foot(s,'給定完整 N，就能決定 f 與 b。IL11 預設關閉，飽和函數為教學假設。');
}
// 20
{
 const s=slide('固定 M 的連乘有適用條件','完整模型的治療與 IL11 使 M 隨隨機狀態改變',
 '若固定每類繁殖分布、突變矩陣與外部環境，就有E[Nt]=N0 M的t次方。完整模型對狀態改變M，因此無條件期望是E[Nt M(St)]，不能一般地將M內的隨機狀態用平均狀態替代。這不妨礙完整P固定，因為P的一列對應該完整狀態，已包含相應環境。區分固定巨大P與隨當前狀態變化的小型M。');
 line(s,'固定環境、固定繁殖規則時',225,75,1130,29);
 eq(s,'E[Nₜ] = N₀ Mᵗ',285,47);
 line(s,'完整模型的一般形式',410,75,1130,29);
 eq(s,'E[Nₜ₊₁] = E[Nₜ M(Sₜ)]',470,43);
 foot(s,'完整 P 仍可固定。不同狀態的 IL11 與治療效果，分別寫在不同的列。');
}
// 21
{
 const s=slide('程式透過抽樣，避免建立巨大 P','一條路徑是一個案例，多次試驗才能描述隨機變異',
 '完整潛在基因型2的22次方為四百一十九萬四千三百零四。但完整群體狀態還須列出這些基因型的各個可能細胞數和治療開關，數量遠大於基因型數。程序只保存實際出現的基因型，並用二項抽樣批次處理計數。默認30次重復，用中位數與10到90百分位表達試驗間變異。所有結果是教學假設下的預測，未經臨床校准。', 'modeling_cancer.py:147-169, 639-760');
 eq(s,'22 位元：2²² = 4,194,304 個潛在基因型',220,35);
 line(s,'population = {基因型: 細胞數}',340,75,1130,33);
 line(s,'simulate() 抽出下一代，不列出整個轉移矩陣',440,75,1130,31);
 line(s,'預設 30 次試驗，呈現中位數與 10–90% 區間',540,75,1130,29);
 foot(s,'百分位區間表示試驗間變異。達計算上限後的後續數值視為缺失。');
}
// Results: freshly reproduced 30-trial experiment groups.
{
 const r=results.therapy;
 const s=slide('結果 1　治療後負荷下降，抗藥族群再擴張','主治療情境：TSG two-hit，IL11 OFF，每組 30 次試驗',
 '資料來自重新執行後的 trajectory.csv，與原有結果逐位元一致。圖顯示第零至三十代，這個時間範圍每代都還有三十次有效觀察。先在原始細胞數上計算中位數及百分位，再做 log10(N+1) 單調轉換。敏感與抗藥曲線是各自中位數，兩者中位數不必相加等於總量中位數。治療各自達到兩萬後開始，起始世代中位數為十五，沒有共同起始線。最低負荷813.5是各試驗nadir的中位數，並非中位曲線最低點。選擇壓力是機制解釋，不能把此模型直接當成病人曲線。',
 'outputs/trajectory.csv, therapy_replicates.csv, run_metadata.json');
 const log=v=>v.map(x=>Math.log10(x+1));
 scientificChart(s,'line',r.generation,[curve('總量 P10',log(r.total_p10),'#A5B4BD',true),
  curve('總量中位',log(r.total_median),INK),curve('總量 P90',log(r.total_p90),'#A5B4BD',true),
  curve('敏感中位',log(r.sensitive_median),BLUE),curve('抗藥中位',log(r.resistant_median),RED)],'log₁₀(細胞數+1)');
 side(s,'813.5','各試驗最低負荷\n的中位數\n\n下降後仍可復發');
 foot(s,'圖 1　前 30 代各代 n=30。實線為中位數，總量虛線為 P10/P90，表示試驗間變異。');
}
{
 const r=results.therapy;
 const s=slide('結果 2　抗藥比例在治療後迅速升高','仍存活且有觀察值的腫瘤，各代抗藥細胞百分比',
 '同樣只顯示前30代，各代有30個活腫瘤觀察值。第十五代中位抗藥比例約百分之一點一六，第十八代約百分之九十點五二，第二十代約百分之九十九點九三。治療先大量移除敏感者，也容許抗藥者持續繁殖。因此比例上升本身不表示治療提高突變率，突變率在程式固定。這是組成選擇以及抗藥數量增長共同結果。若未來出現滅絕，空腫瘤不定義比例，不能把空族群抗藥比例填零。',
 'outputs/trajectory.csv, therapy_trajectories.csv');
 scientificChart(s,'line',r.generation,[curve('P10',r.resistant_percent_p10,'#DBB1A5',true),
  curve('中位數',r.resistant_percent_median,RED),curve('P90',r.resistant_percent_p90,'#DBB1A5',true)],'抗藥細胞 (%)',{max:100});
 side(s,'99.93%','第 20 代\n抗藥比例中位數\n\n治療改變選擇條件');
 foot(s,'圖 2　前 30 代 n=30。虛線為 P10/P90。治療未提高突變率，比例以存活腫瘤為條件。');
}
{
 const a=results.hits['TSG one-hit'], b=results.hits['TSG two-hit'];
 const s=slide('結果 3　one-hit 組較早獲得 TSG 優勢','不施加治療，觀察 18 代，每種門檻 30 次試驗',
 '圖來自one_hit_vs_two_hit.csv。第十八代one-hit中位94617.5，two-hit73879，兩組中位數之比為1.2807，約高28.1%。這是組別中位數比，不是逐細胞或逐試驗因果效應。兩組使用不同種子區段，其他生長與突變參數相同。百分位重疊反映隨機變異，不能推論每一個one-hit試驗都比two-hit大。one-hit與two-hit各自的TSG達標比例使用不同門檻，不能當成相同突變狀態的頻率直接比較。',
 'outputs/one_hit_vs_two_hit.csv, one_hit_vs_two_hit_summary.csv');
 scientificChart(s,'line',a.generation,[curve('one-hit P10',a.total_p10,'#ABC9D9',true),curve('one-hit 中位',a.total_median,BLUE),
 curve('one-hit P90',a.total_p90,'#ABC9D9',true),curve('two-hit P10',b.total_p10,'#DBB1A5',true),
 curve('two-hit 中位',b.total_median,RED),curve('two-hit P90',b.total_p90,'#DBB1A5',true)],'細胞數');
 side(s,'+28.1%','第 18 代\n兩組中位負荷比\n\n94,617.5 / 73,879');
 foot(s,'圖 3　每組 n=30，虛線為 P10/P90。比較功能門檻，其他參數相同，保留隨機重疊。');
}
{
 const a=results.il11['IL11 OFF'], b=results.il11['IL11 ON'];
 const s=slide('結果 4　共享效益提高整體族群負荷','IL11 ON/OFF，不施加治療，每組 30 次試驗',
 '圖由il11_comparison.csv重算百分位。第十八代OFF中位72455.5，ON中位110469.5，兩組中位數比1.52465，約高52.5%。ON生產者比例中位1.603%，共享死亡降低量中位0.026694，約2.67個百分點。ON與OFF維持相同20個基因座及突變機會，關閉時同一位點為passenger。相同seed可重現各次試驗，但族群分歧後不代表逐細胞配對。解釋是本模型的均勻共享效益；缺少分泌成本和空間限制可能放大效益，需作敏感度分析。',
 'outputs/il11_comparison.csv, il11_comparison_summary.csv');
 scientificChart(s,'line',a.generation,[curve('OFF P10',a.total_p10,'#ABC9D9',true),curve('OFF 中位',a.total_median,BLUE),
 curve('OFF P90',a.total_p90,'#ABC9D9',true),curve('ON P10',b.total_p10,'#DBB1A5',true),
 curve('ON 中位',b.total_median,RED),curve('ON P90',b.total_p90,'#DBB1A5',true)],'細胞數',{max:180000});
 side(s,'+52.5%','第 18 代中位比\n\n生產者中位 1.60%\n死亡降低約 2.67\n個百分點');
 foot(s,'圖 4　每組 n=30，虛線為 P10/P90。ON/OFF 保留相同位點與突變機會，機制未校準。');
}
{
 const groups=results.death_rate;
 const s=slide('結果 5　高死亡率組的等大小樣本更多樣','首次跨過 5,000 個細胞後，無放回取樣 5,000 個細胞',
 '每個點代表一次試驗的固定大小樣本Shannon H，以自然對數計算。為了降低跨門檻overshoot影響，從完整跨門檻族群無放回取5000細胞。三組均30/30達標。死亡率0.20、0.35、0.45的H中位分別1.06236、1.80123、2.83497。達標代數中位分別8、14、26.5，因此高死亡率組可能累積更多世代的突變。此比較同時改變周轉與達標時間，只能描述在達標條件下的多樣性關聯。下一步可以同時做固定時間及固定負荷比較，分辨兩者貢獻。橫向抖動僅為避免重疊。',
 'outputs/comparisons.csv; modeling_cancer.py: sample_at_size, run_comparisons');
 const colors=[BLUE,RED,'#66528B'];
 const series=groups.map((g,i)=>({name:`d=${g.death}`,xValues:g.shannon.map((_,j)=>g.death+((j%7)-3)*.002),
  values:g.shannon,fill:colors[i],line:{fill:'none',width:0},marker:{symbol:'circle',size:6}}));
 series.push({name:'中位數',xValues:groups.map(g=>g.death),values:groups.map(g=>g.median_shannon),fill:INK,
  line:{fill:'none',width:0},marker:{symbol:'circle',size:12}});
 const ch=scientificChart(s,'scatter',[],series,'Shannon H',{xTitle:'自然死亡機率 d'});
 ch.xAxis={...ch.xAxis,min:.15,max:.5,majorUnit:.05,numberFormatCode:'0.00'};
 side(s,'1.06\n1.80\n2.83','三組 H 中位數\n\n達標世代中位數\n8、14、26.5');
 foot(s,'圖 5　每組 30/30 達標。彩點為個別試驗，黑色大點為中位數，抖動僅改善重疊。');
}
{
 const r=results.therapy;
 const s=slide('結果 6　這批試驗均觀察到復發','復發：先低於治療起始負荷，再回到該負荷',
 '30/30治療試驗都在觀察期間復發，這只適用於此組設定與種子。復發代數中位28，P10為26，P90為29.1。治療開始代數中位15。所有試驗最後達細胞計算上限，不能以此估計長期生物上限或腫瘤穩態。此模型抗藥細胞治療下的基準期望子代數1.37088仍大於1，為再擴張提供機制支持。不能把30/30解釋為任何真實療法復發概率百分之百。後續可以掃描抗藥殺傷率並記錄滅絕和復發事件。',
 'outputs/therapy_replicates.csv, run_metadata.json');
 scientificChart(s,'bar',r.relapse_histogram.generation,[{name:'復發試驗數',values:r.relapse_histogram.count,fill:RED,valuesFormatCode:'0'}],
  '試驗數',{xTitle:'首次復發世代',labels:true});
 side(s,'30/30','觀察到復發\n\n復發世代中位 28\nP10/P90 = 26/29.1');
 foot(s,'圖 6　n=30。柱高為試驗數。此結論受教學參數與觀察窗口限制，沒有臨床校準。');
}
{
 const s=slide('結果支持的機制與下一步實驗','以可量測的指標檢查解釋，保留模型與臨床的界線',
 '表格把目前觀察、候選機制與後續實驗連在一起。抗藥殺傷率由0.04掃描至更高值，固定其他設定並做多次試驗，記錄復發比例、事件時間及滅絕。對IL11掃描共享效益上限與半飽和參數，加入分泌成本或空間限制後檢查52.5%的差異是否仍存在。對多樣性同時固定世代與固定負荷，並呈現未達標試驗。每項假說都需要獨立的控制實驗。未經校准不能推出最佳臨床治療策略。');
 table(s,[['觀察','可能機制','下一步與指標'],['抗藥佔優與再擴張','敏感與抗藥的存活差異','掃描抗藥殺傷率\n復發比例、事件時間'],
  ['IL11 ON 負荷較高','均勻共享死亡降低','掃描 b_max、K 與分泌成本\n負荷差異、生產者比例'],
  ['高死亡組 H 較高','周轉與達標時間共同影響','固定世代及負荷分別比較\nH、達標率與代數']],70,215,1140,355,[300,350,490],25);
 foot(s,'各圖描述指定參數下的結果。後續需做敏感度分析與資料校準，再評估生物適用性。');
}
// Four-formula synthesis, then mathematical appendices.
{
 const s=slide('四個公式的完整讀法','狀態、抽樣、機率分布與條件期望',
 '回到一開始四個公式。第一行定義完整狀態。第二行講一次隨機模擬如何產生下一步。第三行描述所有可能狀態的機率分布。第四行只計算下一代細胞數期望。聽眾應能解釋為什麼pi在P左邊，為什麼P和M的索引不同，以及為何治療開關不能省略。接下來附錄給完整數學形式。');
 const rows=[['Sₜ = (Nₜ,Aₜ)','現在有哪些細胞？治療是否已啟動？'],['Sₜ₊₁ | Sₜ=s ∼ P(s,·)','給定現在，抽取下一個完整狀態'],['πₜ = π₀ Pᵗ','加總所有可能路徑，得到狀態分布'],['E[Nₜ₊₁ | Sₜ=s] = nM(s)','給定現在，計算下一代平均細胞數']];
 rows.forEach(([a,b],i)=>{eq(s,a,225+i*92,35,75,550,70);text(s,b,655,236+i*92,550,70,27);});
}
// 23 appendix
{
 const s=slide('附錄 A　完整轉移的生長參數','以下公式使用預設係數，並在停止條件尚未觸發時適用',
 '定義IO為oncogene達門檻的指標，IT為任一TSG功能失活達門檻的指標。每類優勢最多給一次。d取零下限。alpha是這次轉移使用的治療開關，與門檻檢查前的a區分。自然生存與治療生存相乘得到有效死亡delta。治療功能關閉時alpha等於a且初始化a零。空群體IL11 benefit零。', 'modeling_cancer.py:123-144, 230-244');
 eq(s,'d_g(N) = max{0, 0.30−0.08 I_O(g)−0.08 I_T(g)−b(N)}',215,31);
 eq(s,'q_g = 0.04+0.12 I_O(g)+0.12 I_T(g)',310,35);
 eq(s,'α = a ∨ 1{治療功能開啟且 ∑g n_g ≥ 20,000}',405,33);
 eq(s,'δ_g(s) = 1−(1−d_g(N))(1−α k_g)',500,36);
 foot(s,'I_O、I_T 是是否達標的指標。敏感 k_g = 0.90，抗藥 k_g = 0.04。');
}
// 24 appendix
{
 const s=slide('附錄 B　獨立位元突變與子代分配','L = 22，μ = 0.001，突變不可逆',
 'g和h視為已突變位元集合，|g|是突變數。如果g是h子集，那麼h新增的位元數為|h|-|g|，其餘仍未突變位元數為L-|h|，乘上獨立概率即得到Qgh。g不是h子集意味著回復突變，在此模型概率為零。給定Bg，所有來自g父細胞的子代按Q多項分布。最後對所有父基因型求和。', 'modeling_cancer.py:147-169');
 eq(s,'Q_g,h = μ^(|h|−|g|) (1−μ)^(L−|h|)，若 g ⊆ h',215,33);
 line(s,'其餘情況 Q_g,h = 0，每個 row 的總和為 1',310,75,1130,29);
 eq(s,'(Y_g,h)ₕ | B_g ∼ Multinomial(B_g, (Q_g,h)ₕ)',410,34);
 eq(s,'Nₕ,ₜ₊₁ = ∑g Y_g,h',510,43);
 foot(s,'Y_g,h 是來自父基因型 g、最後成為子基因型 h 的細胞數。');
}
// 25 appendix
{
 const s=slide('附錄 C　用生成函數定義精確 P','取出指定次方的係數，即得到指定下一代族群的機率',
 'u_g是單個子代基因型的生成函數。Fg是單個父細胞全部後代的多變量生成函數，分別包含死亡、二子代與三子代。給定狀態父細胞獨立，因此全族群函數H是各Fg的n_g次方乘積。Pss′是H中對應n′向量的系數，再乘治療開關一致的指標。方括號[z^n′]表示系數提取。投影片為方便閱讀把多變量符號z和n寫成簡潔形式，它們都是向量，z^n′表示所有h的z_h^n′_h乘積。');
 eq(s,'u_g(z) = ∑h Q_g,h z_h',210,36);
 eq(s,'F_g(z|s) = δ_g+(1−δ_g)[(1−q_g)u_g²+q_g u_g³]',300,32);
 eq(s,'H_s(z) = ∏g F_g(z|s)^n_g',405,37);
 eq(s,'P_s,s′ = 1{a′=α} [z^n′] H_s(z)',505,37);
 foot(s,'z、n′ 為向量。[z^n′] 代表取出 ∏h z_h^n′_h 的係數。');
}
// 26 appendix
{
 const s=slide('附錄 D　路徑機率與停止時間','滅絕是生物吸收，計算上限則是觀察停止',
 '單條路徑概率等於初始狀態概率乘沿途轉移概率。tau是第一次出現停止條件的世代。最後觀測代數可通過把t加入狀態形成時間齊次的執行鏈。若要逐步識別復發，需額外保存治療開始負荷以及曾經下降的標記，否則當前N和A不足以判定該歷史事件。代碼事後從history判斷復發，不在simulate內切換模式。', 'modeling_cancer.py:218-229, 258-279');
 eq(s,'Pr(s₀,…,s_T) = π₀(s₀) ∏ₜ P_sₜ,sₜ₊₁',215,35);
 eq(s,'τ = inf{t : 滅絕、取樣達標、計算上限或 t=G}',320,32);
 line(s,'含停止控制的狀態可擴充為 Xₜ = (Nₜ,Aₜ,t,Cₜ)',435,75,1130,29);
 line(s,'逐步判讀復發，還需保存治療起始負荷與曾下降標記',535,75,1130,28);
 foot(s,'Cₜ 記錄執行中或停止原因。未觀察到復發不等於永久治癒。');
}

const candidate=path.join(build,'candidate.pptx');
await (await PresentationFile.exportPptx(deck)).save(candidate);
await fs.writeFile(path.join(build,'slides.json'),JSON.stringify(records,null,2));
const finalPath=path.join(out,process.env.PRESENTATION_FILENAME??'tumor_markov_results_zh_tw_v3.pptx');
const result=await finalizePresentation({workspaceDir:workspace,candidatePath:candidate,finalPath,
 pythonExecutable:process.env.RUNTIME_PYTHON??path.join(runtime,'python/bin/python3'),
 integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit',
   ...[...new Set(nativeTables)].flatMap(n=>['--require-native-table-slide',String(n)])],
 explicitTotalSlideCount:records.length,requiredNativeTableOwnerSlides:[...new Set(nativeTables)],
 requiredNativeChartOwnerSlides:chartOwners,materializeLiteralChartWorkbooks:true,nativeChartTargetApplication:'powerpoint',
 fontPolicy:{basis:'design',families:[FONT]},verifyArtifactToolImport:true,
 receiptPath:path.join(build,`${path.basename(finalPath)}.${Date.now()}.validation.json`)});
console.log(JSON.stringify({finalPath,result}));

// Render the checked final file, rather than using draft previews as proof.
const checked=await PresentationFile.importPptx(await FileBlob.load(finalPath));
for(let i=0;i<records.length;i++){
 const s=checked.slides.items[i];
 const blob=await checked.export({slide:s,format:'png',scale:1.25});
 await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await blob.arrayBuffer()));
 const layout=await s.export({format:'layout'});
 await fs.writeFile(path.join(build,`slide-${String(i+1).padStart(2,'0')}.layout.json`),await layout.text());
}
console.log(`Rendered all ${records.length} final slides.`);
