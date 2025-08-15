from gpt4 import save
from openai import OpenAI
import datetime, inspect, os, json, random, sys, codecs
sys.stdout.reconfigure(encoding='utf-8')

MODELS = [
    "o4-mini",
    "o4-mini-deep-research",
    "o3-pro",
    "o3",
    "o3-deep-research",
    "o1-pro",
    "o1",
    "gpt-4.1",
    "gpt-4.1-mini",
    "gpt-4.1-nano",
    "gpt-4o",
    "gpt-4o-mini",
]

MODEL = MODELS[0] # CHANGE INDEX FOR DIFFERENT MODEL

SCREENSHOTS = [
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543577957044264/00-Sphere.png?ex=68821b62&is=6880c9e2&hm=fe8430de6d69c5b391712558e6dc64189659b44ef9f58d5d374869c58f305478&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543578213027900/01-Cube.png?ex=68821b62&is=6880c9e2&hm=a0155828b41149b34b8b017c6a1beae6b61bf28656a8143db481a76088d9e1cb&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543578477002753/02-Square.png?ex=68821b62&is=6880c9e2&hm=9eb287a78e2dfe343ce725d7dd7b9cd8a60cb78cad9fcaca973ab6afbd14ef05&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543578732990516/03-Cylinder.png?ex=68821b62&is=6880c9e2&hm=8ef91923c2b47c65c361215b2ef093698d698ef0cc1a6946f82b5b3d56f60855&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543579022528522/04-Two_Spheres.png?ex=68821b62&is=6880c9e2&hm=05dbc6b3c0751f478cc6fc6ec293db7dd2e2d022dadd72d0e8b9a4e9038672c8&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543579274182717/05-Disc.png?ex=68821b62&is=6880c9e2&hm=0c2c95490eb28b3f9ba0e509510089676aa54973c0f285b91d9095f611595f5a&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543579886424075/06-Quad_with_Circle_Hole.png?ex=68821b62&is=6880c9e2&hm=eaa8d2624e2bd25330ba727d4878679d82eca538cd752b426bdb5f8c88775837&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543580129824860/07-Plane_with_Triangle_Hole.png?ex=68821b62&is=6880c9e2&hm=79d4d17be01f8ab47752aa29dc3544eab8c125d42de544649515a2ecf1399d05&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543580373090464/08-Ellipsoid.png?ex=68821b62&is=6880c9e2&hm=b927f122f262b41c864fece3e1f1a030039647645ef39fd3b29f81f13dcdbb47&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543577676152853/09-Cut_Sphere.png?ex=68821b62&is=6880c9e2&hm=0c17d6853d2260224ee548c7aa2a99321a8055ab8cd5d63937494d42f7ee97e6&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543599893250048/10-Sphere_with_Sphere_Hole.png?ex=68821b67&is=6880c9e7&hm=ba7c75d24cd70bbdefe57025ca6c67aa7e33ebe16cbe4d90592ae276cf2038bc&",
    "https://cdn.discordapp.com/attachments/1203556084133138473/1397543599544991888/11-Sphere_with_Cylinder_Hole.png?ex=68821b67&is=6880c9e7&hm=d8b0b0dc4248c4ce13fc02e7dc4ed915be9ab2f5e24d4838024671650df69351&",
]

file = open(f"dump_{MODEL}_screenshots.txt", "w")
client = OpenAI()
res = []
NUMBER_OF_RUNS = 10

for k in range(len(SCREENSHOTS)):
    res.append([])
    for i in range(NUMBER_OF_RUNS):
        content = [
            {"type": "input_text", "text": "What is the solution to this question?"},
            {
                "type": "input_image",
                "image_url": SCREENSHOTS[k],
            },
        ]
        try:
            response = client.responses.create(
                model = MODEL,
                input = [{
                    "role": "user",
                    "content": content
                }],
            ).output_text.strip().replace("\\n", "\n")
            content.append({"role": "assistant", "content": response})
            now = datetime.datetime.now().strftime("%Y.%m.%d-%H.%M.%S.%f")
            caller_location = os.path.dirname(os.path.abspath(inspect.stack()[0].filename))
            if not os.path.exists(f"{caller_location}/GPT-4 History"):
                os.mkdir(f"{caller_location}/GPT-4 History")
            file_history = open(f"{caller_location}/GPT-4 History/{now}[{int(round(random.random(), 6) * 1000000)}]Screenshot_Shape{k}_Run{i}.log", "w")
            json.dump(content, file_history, indent=4)
            file_history.close()
            print(f"{k}.{i}")
        except Exception as e:
            response = "ERROR"
            print(e)
        res[-1].append(response)

json.dump(res, file, indent=4)
file.close()
