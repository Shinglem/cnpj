package io.github.shinglem.cnpj.webserver


import io.netty.handler.codec.http.HttpHeaderNames
import io.vertx.core.http.HttpMethod
import io.vertx.core.json.JsonObject
import io.vertx.ext.web.Router
import io.vertx.ext.web.handler.BodyHandler
import io.vertx.ext.web.handler.LoggerHandler
import io.vertx.ext.web.handler.StaticHandler
import io.vertx.kotlin.core.http.listenAwait
import io.vertx.kotlin.coroutines.CoroutineVerticle
import kotlinx.coroutines.launch
import org.slf4j.LoggerFactory


class WebserverVerticle() : CoroutineVerticle() {


    private val logger = LoggerFactory.getLogger(this::class.java.name)


    var port: Int = 12000


    override suspend fun start() {
        try {
            logger.debug("---------web start---------" + this.deploymentID)

            val server = vertx.createHttpServer()

            val router = Router.router(vertx)




            router.route().order(-1).handler(LoggerHandler.create())
            router.get("/web/*").order(-1).handler(StaticHandler.create("/web").setCachingEnabled(false))

            router.route("/api/*").order(-1).handler(BodyHandler.create())



            router.errorHandler(404) {
                logger.info("[404 handler]")
                it.response().statusCode = 404
                if (it.request().method() != HttpMethod.HEAD) {
                    // If it's a 404 let's send a body too
                    it.response()
                        .putHeader(HttpHeaderNames.CONTENT_TYPE, "application/json; charset=\"utf-8\"")
                        .end(responseBuild(ResponseCode.RESULT_NOK, "404 not found").encodePrettily())
                } else {
                    it.response().end()
                }
            }

            router.errorHandler(500) {
                logger.info("[500 handler]");
                it.response().statusCode = 500
                if (it.request().method() != HttpMethod.HEAD) {
                    // If it's a 404 let's send a body too
                    it.response()
                        .putHeader(HttpHeaderNames.CONTENT_TYPE, "application/json; charset=\"utf-8\"")
                        .end(responseBuild(ResponseCode.RESULT_NOK, "Internal Server Error").encodePrettily())
                } else {
                    it.response().end()
                }
            }

            router.errorHandler(403) {
                logger.info("[403 handler]");
                it.response().statusCode = 403
                if (it.request().method() != HttpMethod.HEAD) {
                    // If it's a 404 let's send a body too
                    it.response()
                        .putHeader(HttpHeaderNames.CONTENT_TYPE, "application/json; charset=\"utf-8\"")
                        .end(responseBuild(ResponseCode.RESULT_NOK, "403 forbiden").encodePrettily())
                } else {
                    it.response().end()
                }
            }

            router.route().order(0).handler { rc ->
                logger.info("[rootRouter]")

                val response = rc.response()
                response.isChunked = true
                response.putHeader("content-type", "application/json; charset=\"utf-8\"")
                rc.next()
            }


            router.route().failureHandler {
                logger.info("[failure Handler]");


                val response = it.response()
                response.putHeader("content-type", "application/json; charset=\"utf-8\"")

                val errorCode = it.statusCode()



                when (errorCode) {

                    404 -> {
                        logger.info("[404 fail]");
                        it.response().statusCode = 404
                        if (it.request().method() != HttpMethod.HEAD) {
                            // If it's a 404 let's send a body too
                            it.response()
                                .putHeader(HttpHeaderNames.CONTENT_TYPE, "application/json; charset=\"utf-8\"")
                                .end(responseBuild(ResponseCode.RESULT_NOK, "404 not found").encodePrettily())
                        } else {
                            it.response().end()
                        }
                    }
                    500 -> {
                        logger.info("[500 fail]");
                        it.response().statusCode = 500
                        if (it.request().method() != HttpMethod.HEAD) {
                            // If it's a 404 let's send a body too
                            it.response()
                                .putHeader(HttpHeaderNames.CONTENT_TYPE, "application/json; charset=\"utf-8\"")
                                .end(responseBuild(ResponseCode.RESULT_NOK, "Internal Server Error").encodePrettily())
                        } else {
                            it.response().end()
                        }
                    }
                    403 -> {
                        logger.info("[403 fail]");
                        it.response().statusCode = 403
                        if (it.request().method() != HttpMethod.HEAD) {
                            // If it's a 404 let's send a body too
                            it.response()
                                .putHeader(HttpHeaderNames.CONTENT_TYPE, "application/json; charset=\"utf-8\"")
                                .end(responseBuild(ResponseCode.RESULT_NOK, "403 forbiden").encodePrettily())
                        } else {
                            it.response().end()
                        }
                    }
                    else -> {
                        logger.info("[$errorCode fail]");
                        it.response().statusCode = 500
                        if (it.request().method() != HttpMethod.HEAD) {
                            // If it's a 404 let's send a body too
                            it.response()
                                .putHeader(HttpHeaderNames.CONTENT_TYPE, "application/json; charset=\"utf-8\"")
                                .end(responseBuild(ResponseCode.RESULT_NOK, "Internal Server Error").encodePrettily())
                        } else {
                            it.response().end()
                        }
                    }
                }


            }





            server
                .requestHandler(router)
                .listenAwait(port)
            logger.info("web start success , listening ${port}")

            val httpClient = vertx.createHttpClient()
            httpClient.getAbs("http://127.0.0.1:12000/web/index.html"){
                logger.info("get web success ")
            }.end()

        } catch (e: Exception) {

            logger.error("web start fail ", e)
        }
    }

}

fun responseBuild(code: ResponseCode, msg: String): JsonObject = responseBuild(code, msg, JsonObject())

fun responseBuild(code: ResponseCode, msg: String, data: JsonObject): JsonObject {
    return JsonObject().put("code", code.code).put("msg", msg).put("data", data)
}

enum class ResponseCode(val code: Int) {
    RESULT_OK(0), RESULT_NOK(999)
}